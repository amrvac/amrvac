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

    integer :: igrid,iigrid

    call init_comm_fix_conserve(1,ndim,^ND)

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ! Make zero the magnetic fluxes
       ! Fake advance, storing electric fields at edges
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
       call faces2centers(ixG^LL,ps(igrid))
    end do
    !$OMP END PARALLEL DO
  end subroutine recalculateB

  !> fake advance a step to calculate magnetic field
  subroutine fake_advance(igrid,idim^LIM,s)
    use mod_global_parameters
    use mod_fix_conserve

    integer       :: igrid,idim^LIM
    type(state)   :: s

    double precision             :: dx^D
    double precision             :: fC(ixG^T,1:nwflux,1:ndim)
    double precision             :: fE(ixG^T,1:ndir)

    dx^D=rnode(rpdx^D_,igrid);
    !call set_tmpGlobals(igrid)

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
    double precision             :: fE(ixI^S,1:ndir)
    double precision             :: dx^D

    integer                            :: ixIs^L,ixO^L,idir
    double precision                   :: xC(ixI^S,1:ndim), A(ixI^S,1:ndir)
    double precision                   :: circ(ixI^S,1:ndim), dxidir

    associate(ws=>s%ws,x=>s%x)

    A(:^D&,:)=zero
    ws(:^D&,:)=zero

    ixIs^L=s%ixGs^L;
    ixO^L=ixI^L^LSUBnghostcells;

    call b_from_vectorpotentialA(ixIs^L, ixI^L, ixO^L, ws, x, A)

    ! This is important only in 3D
    do idir=1,ndim
       fE(ixI^S,idir) =-A(ixI^S,idir)*dxlevel(idir)
    end do

    end associate
  end subroutine fake_update

  !> calculate cell-center values from face-center values
  subroutine faces2centers(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s
    integer                            :: hxO^L, idim

    associate(w=>s%w, ws=>s%ws)

    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxO^L=ixO^L-kr(idim,^D);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixO^S,iw_mag(idim))=half/s%surface(ixO^S,idim)*&
        (ws(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
        +ws(hxO^S,idim)*s%surfaceC(hxO^S,idim))
    end do

    end associate

  end subroutine faces2centers

  !> calculate cell-center values from face-center values in 4th order
  subroutine faces2centers4(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer                            :: gxO^L, hxO^L, jxO^L, idim

    associate(w=>s%w, ws=>s%ws)

    do idim=1,ndim
      gxO^L=ixO^L-2*kr(idim,^D);
      hxO^L=ixO^L-kr(idim,^D);
      jxO^L=ixO^L+kr(idim,^D);

      ! Interpolate to cell barycentre using fourth order central formula
      w(ixO^S,iw_mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
             ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
         +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
         +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
               -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    end do

    end associate

  end subroutine faces2centers4

  !> calculate cell-center values from face-center values in 6th order
  subroutine faces2centers6(ixO^L,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixO^L
    type(state)                        :: s

    integer                            :: fxO^L, gxO^L, hxO^L, jxO^L, kxO^L, idim

    associate(w=>s%w, ws=>s%ws)

    do idim=1,ndim
      fxO^L=ixO^L-3*kr(idim,^D);
      gxO^L=ixO^L-2*kr(idim,^D);
      hxO^L=ixO^L-kr(idim,^D);
      jxO^L=ixO^L+kr(idim,^D);
      kxO^L=ixO^L+2*kr(idim,^D);

      ! Interpolate to cell barycentre using sixth order central formula
      w(ixO^S,iw_mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
         (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
           -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
          +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
          +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
           -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
            +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    end do

    end associate

  end subroutine faces2centers6

  !> calculating the divergence of a face-allocated vector field.
  subroutine div_staggered(ixO^L,s,divv)
    use mod_global_parameters

    integer, intent(in)           :: ixO^L
    type(state)                   :: s
    double precision              :: divv(ixO^S)

    integer                       :: hxO^L,idim

    divv=zero

    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
       hxO^L=ixO^L-kr(idim,^D);
       ! Calculate divergence by taking differences of fluxes
       divv(ixO^S)=divv(ixO^S)+(s%ws(ixO^S,idim)*s%surfaceC(ixO^S,idim)&
                               -s%ws(hxO^S,idim)*s%surfaceC(hxO^S,idim))&
                               /s%dvolume(ixO^S)
    end do

  end subroutine div_staggered

  !> calculate magnetic field from vector potential
  subroutine b_from_vectorpotential(ixIs^L, ixI^L, ixO^L, ws, x)
    use mod_global_parameters

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws)
    double precision, intent(in)       :: x(ixI^S,1:ndim)

    double precision                   :: Adummy(ixI^S,1:ndir)

    call b_from_vectorpotentialA(ixIs^L, ixI^L, ixO^L, ws, x, Adummy)

  end subroutine b_from_vectorpotential

  !> calculate magnetic field from vector potential A at cell edges
  subroutine b_from_vectorpotentialA(ixIs^L, ixI^L, ixO^L, ws, x, A)
    use mod_global_parameters
    use mod_usr_methods, only: usr_init_vector_potential

    integer, intent(in)                :: ixIs^L, ixI^L, ixO^L
    double precision, intent(inout)    :: ws(ixIs^S,1:nws),A(ixI^S,1:ndir)
    double precision, intent(in)       :: x(ixI^S,1:ndim)

    integer                            :: ixC^L, hxC^L, ixCp^L, ixCm^L, hxO^L, idim, idim1, idim2, idir
    double precision                   :: xC(ixI^S,1:ndim)
    double precision                   :: circ(ixI^S,1:ndim)

    A=zero
    ws=zero

    {ixCmax^D=ixOmax^D;}
    {ixCmin^D=ixOmin^D-1;} ! Extend range by one
    do idir=7-2*ndim,ndir
      do idim=1,ndim
        ! Get edge coordinates
        if (idim/=idir) then
          xC(ixC^S,idim)=x(ixC^S,idim)+half*block%dx(ixC^S,idim)
        else
          xC(ixC^S,idim)=x(ixC^S,idim)
        end if
      end do
      ! Initialise vector potential at the edge
      call usr_init_vector_potential(ixI^L, ixC^L, xC, A(ixI^S,idir), idir)
    end do

    ! Set NaN to zero (can happen e.g. on axis):
    where(A(ixI^S,1:ndir)/=A(ixI^S,1:ndir))
       A(ixI^S,1:ndir)=zero
    end where

    ! sub integrals A ds
    A(ixC^S,1:ndir)=A(ixC^S,1:ndir)*block%dsC(ixC^S,1:ndir)

    ! Take the curl of the vector potential 
    circ(:^D&,:) = zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=1,ndir ! Direction of line integral
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

  end subroutine b_from_vectorpotentialA

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
       
!    case (limiter_weno5)
!       call WENO5limitervar(ixI^L,ixC^L,idir,q,qL,qR)
!       
!    case (limiter_wenoZP)
!       call WENOZPlimitervar(ixI^L,ixC^L,idir,dxlevel(idir),q,qL,qR)
       
    case default

       dqC(gxC^S)= qR(gxC^S)-qL(gxC^S)
       call dwlimiter2(dqC,ixI^L,gxC^L,idir,typelimiter,ldq,rdq)
       qL(ixC^S) = qL(ixC^S) + half*ldq(ixC^S)
       qR(ixC^S) = qR(ixC^S) - half*rdq(jxC^S)
    end select

  end subroutine reconstruct

end module mod_constrained_transport
