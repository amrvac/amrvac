module mod_constrained_transport
  implicit none
  public

contains
  !> re-calculate the magnetic field from the vector potential in a completely
  !> divergency free way
  subroutine recalculateB
    use mod_global_parameters
    use mod_fix_conserve
    use mod_physics
    use mod_ghostcells_update

    integer :: igrid,iigrid

    call init_comm_fix_conserve(1,ndim,3)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! Make zero the magnetic fluxes
       ! Fake advance, storing electric fields at edges
       block=>ps(igrid)
       call fake_advance(igrid,1,3,ps(igrid))

    end do
    !$OMP END PARALLEL DO

    ! Do correction
    call recvflux(1,ndim)
    call sendflux(1,ndim)
    call fix_conserve(ps,1,ndim,1,3)

    call fix_edges(ps,1,3)

    ! Now we fill the centers for the staggered variables
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
          ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
       ! update cell center magnetic field
       call phys_face_to_center(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
          ps(igrid))
       call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
          ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
    end do
    !$OMP END PARALLEL DO

    call getbc(global_time,0.d0,ps,iwstart,nwgc)

  end subroutine recalculateB

  !> fake advance a step to calculate magnetic field
  subroutine fake_advance(igrid,idimmin,idimmax,s)
    use mod_global_parameters
    use mod_fix_conserve

    integer       :: igrid,idimmin,idimmax
    type(state)   :: s

    double precision             :: dx1,dx2,dx3
    double precision             :: fC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:nwflux,1:ndim)
    double precision             :: fE(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,sdim:3)

    dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);

    call fake_update(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,s,fC,fE,dx1,dx2,&
       dx3)

    call store_flux(igrid,fC,idimmin,idimmax,3)
    call store_edge(igrid,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,fE,idimmin,&
       idimmax) 

  end subroutine fake_advance

  !>fake update magnetic field from vector potential
  subroutine fake_update(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,s,fC,&
     fE,dx1,dx2,dx3)
    use mod_global_parameters

    integer       :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    type(state)   :: s
    double precision             :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision             :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)
    double precision             :: dx1,dx2,dx3

    integer                            :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
       ixIsmax2,ixIsmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir
    double precision                   :: A(s%ixGsmin1:s%ixGsmax1,&
       s%ixGsmin2:s%ixGsmax2,s%ixGsmin3:s%ixGsmax3,1:3)

    associate(ws=>s%ws,x=>s%x)

    A=zero
    ws=zero

    ixIsmin1=s%ixGsmin1;ixIsmin2=s%ixGsmin2;ixIsmin3=s%ixGsmin3
    ixIsmax1=s%ixGsmax1;ixIsmax2=s%ixGsmax2;ixIsmax3=s%ixGsmax3;
    ixOmin1=ixImin1+nghostcells;ixOmin2=ixImin2+nghostcells
    ixOmin3=ixImin3+nghostcells;ixOmax1=ixImax1-nghostcells
    ixOmax2=ixImax2-nghostcells;ixOmax3=ixImax3-nghostcells;

    fC=0.d0
    call b_from_vector_potentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,&
       ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, A)

    ! This is important only in 3D
    do idir=sdim,3
       fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          idir) =-A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir)
    end do

    end associate
  end subroutine fake_update

  !> calculate magnetic field from vector potential A at cell edges
  subroutine b_from_vector_potentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
     ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, A)
    use mod_global_parameters
    use mod_usr_methods, only: usr_init_vector_potential

    integer, intent(in)                :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
       ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:nws),A(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:3)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer                            :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
       ixCmax2,ixCmax3, hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3, idim1,&
        idim2, idir
    double precision                   :: xC(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:ndim),xCC(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:ndim)
    double precision                   :: circ(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:ndim)

    A=zero
    ! extend one layer of cell center locations in xCC
    xC=0.d0
    xCC=0.d0
    xCC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)=x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    
    xCC(ixIsmin1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=x(ixImin1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    xCC(ixIsmin1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1)=x(ixImin1,ixImin2,&
       ixImin3,1)-block%dx(ixImin1,ixImin2,ixImin3,1)
    
    
    xCC(ixImin1:ixImax1,ixIsmin2,ixImin3:ixImax3,1:ndim)=x(ixImin1:ixImax1,&
       ixImin2,ixImin3:ixImax3,1:ndim)
    xCC(ixIsmin1:ixIsmax1,ixIsmin2,ixIsmin3:ixIsmax3,2)=x(ixImin1,ixImin2,&
       ixImin3,2)-block%dx(ixImin1,ixImin2,ixImin3,2)
    
    
    xCC(ixImin1:ixImax1,ixImin2:ixImax2,ixIsmin3,1:ndim)=x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3,1:ndim)
    xCC(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3,3)=x(ixImin1,ixImin2,&
       ixImin3,3)-block%dx(ixImin1,ixImin2,ixImin3,3)
    
    
    xCC(ixImin1:ixImax1,ixIsmin2,ixIsmin3,1)=x(ixImin1:ixImax1,ixImin2,ixImin3,&
       1)
    xCC(ixIsmin1,ixImin2:ixImax2,ixIsmin3,2)=x(ixImin1,ixImin2:ixImax2,ixImin3,&
       2)
    xCC(ixIsmin1,ixIsmin2,ixImin3:ixImax3,3)=x(ixImin1,ixImin2,ixImin3:ixImax3,&
       3)
   

    do idir=sdim,3
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2)
      ixCmin3=ixOmin3-1+kr(idir,3);
      do idim1=1,ndim
        ! Get edge coordinates
        if (idim1/=idir) then
          xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=xCC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)+half*block%dx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idim1)
        else
          xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=xCC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
        end if
      end do
      ! Initialise vector potential at the edge
      call usr_init_vector_potential(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
         ixIsmax2,ixIsmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
          xC, A(ixIsmin1:ixIsmax1,ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,idir),&
          idir)
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)*block%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
    end do

    ! Take the curl of the vector potential 
    circ=zero
    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          if(lvc(idim1,idim2,idir)==0) cycle
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
          hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)+lvc(idim1,idim2,idir)* (A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idir)-A(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
             hxCmin3:hxCmax3,idir))
        end do
      end do
      ! Divide by the area of the face to get B
      where(block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)==0)
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,idim1)
      end where
      ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) = circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

  end subroutine b_from_vector_potentialA

  !> Reconstruct scalar q within ixO^L to 1/2 dx in direction idir
  !> Return both left and right reconstructed values 
  subroutine reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,q,qL,qR)
    use mod_limiter
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
    double precision, intent(in)       :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out)      :: qL(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), qR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    double precision                   :: qC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision,dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: dqC,ldq,rdq
    integer                            :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,gxCmin1,&
       gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3

    jxCmin1=ixCmin1+kr(idir,1);jxCmin2=ixCmin2+kr(idir,2)
    jxCmin3=ixCmin3+kr(idir,3);jxCmax1=ixCmax1+kr(idir,1)
    jxCmax2=ixCmax2+kr(idir,2);jxCmax3=ixCmax3+kr(idir,3);
    gxCmin1=ixCmin1-kr(idir,1);gxCmin2=ixCmin2-kr(idir,2)
    gxCmin3=ixCmin3-kr(idir,3);gxCmax1=jxCmax1;gxCmax2=jxCmax2
    gxCmax3=jxCmax3;
    hxCmin1=gxCmin1+kr(idir,1);hxCmin2=gxCmin2+kr(idir,2)
    hxCmin3=gxCmin3+kr(idir,3);hxCmax1=gxCmax1+kr(idir,1)
    hxCmax2=gxCmax2+kr(idir,2);hxCmax3=gxCmax3+kr(idir,3);

    qR(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = q(hxCmin1:hxCmax1,&
       hxCmin2:hxCmax2,hxCmin3:hxCmax3)
    qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3) = q(gxCmin1:gxCmax1,&
       gxCmin2:gxCmax2,gxCmin3:gxCmax3)

    dqC(gxCmin1:gxCmax1,gxCmin2:gxCmax2,gxCmin3:gxCmax3)= qR(gxCmin1:gxCmax1,&
       gxCmin2:gxCmax2,gxCmin3:gxCmax3)-qL(gxCmin1:gxCmax1,gxCmin2:gxCmax2,&
       gxCmin3:gxCmax3)
    call dwlimiter2(dqC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       gxCmin1,gxCmin2,gxCmin3,gxCmax1,gxCmax2,gxCmax3,idir,&
       type_limiter(block%level),ldq,rdq)
    qL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qL(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3) + half*ldq(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3)
    qR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = qR(ixCmin1:ixCmax1,&
       ixCmin2:ixCmax2,ixCmin3:ixCmax3) - half*rdq(jxCmin1:jxCmax1,&
       jxCmin2:jxCmax2,jxCmin3:jxCmax3)

  end subroutine reconstruct

end module mod_constrained_transport
