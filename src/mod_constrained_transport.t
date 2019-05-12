module mod_constrained_transport
  implicit none
  public

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
      w(ixO^S,iw_mag(idim))=(half*s%dx(ixO^S,idim)/s%dvolume(ixO^S))*&
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

  !> show on screen the divergence of a face-allocated vector field.
  subroutine show_div_staggered(ixO^L,s)
    use mod_global_parameters

    integer, intent(in)           :: ixO^L
    type(state)                   :: s

    double precision              :: divv(ixG^T)
    integer                       :: hxO^L,idim

    print *, 'In div_stg'
    print *, 'x min = ', s%x(ixMlo^D,1:ndim)

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

    call printarray(ixG^LL,ixO^L,divv)

  end subroutine show_div_staggered

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

  !> update faces
  subroutine updatefaces(ixI^L,ixO^L,qdt,fC,fE,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qdt
    type(state)                        :: s
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,1:ndir)

    integer                            :: ix^L(1:ndim)
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2,i,j,k
    double precision                   :: circ(ixI^S,1:ndim)

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax^D=ixOmax^D;
    ixCmin^D=ixOmin^D-1;

    fE(ixI^S,1:ndir)=zero

    do idim1=1,ndim 
      iwdim1 = iw_mag(idim1)
      do idim2=1,ndim
        iwdim2 = iw_mag(idim2)
        do idir=7-2*ndim,ndir ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxC^L=ixC^L+kr(idim1,^D);
            ixCp^L=ixC^L+kr(idim2,^D);
            ! Interpolate to edges
            fE(ixC^S,idir)=qdt*quarter*((fC(ixC^S,iwdim1,idim2)&
            +fC(jxC^S,iwdim1,idim2))/dxlevel(idim1)&
            -(fC(ixC^S,iwdim2,idim1)+fC(ixCp^S,iwdim2,idim1))&
            /dxlevel(idim2))

            if (.not.slab) then
              where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixC^S,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      do idim2=1,ndim
        do idir=1,ndir ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Decrease bottom limit by one
    do idim1=1, ndim
      ixmax^D(idim1)=ixOmax^D;ixmin^D=ixOmin^D-1;!-kr(^D,idim1);
    end do

    ! Divide by the area of the face to get dB/dt

    do idim1=1,ndim
      ixC^L=ix^L(idim1);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate

  end subroutine updatefaces

  !> update faces UCT2
  subroutine update_faces_uct2(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fE,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qdt
    type(state)                        :: s
    double precision, intent(in)       :: vbarC(ixI^S,1:ndir,2)
    double precision, intent(in)       :: cbarmin(ixI^S,ndim)
    double precision, intent(in)       :: cbarmax(ixI^S,ndim)
    double precision, intent(inout)    :: fE(ixI^S,1:ndir)

    double precision                   :: vtilL(ixI^S,2)
    double precision                   :: vtilR(ixI^S,2)
    double precision                   :: btilL(s%ixGs^S,ndim)
    double precision                   :: btilR(s%ixGs^S,ndim)
    double precision                   :: sqrtg(s%ixGs^S)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    integer                            :: ixGs^L
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir
    double precision                   :: circ(ixI^S,1:ndim), dxidir
    integer                            :: i^D

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ixGs^L=s%ixGs^L;

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    fE(ixI^S,1:ndir)=zero

    do idir=7-2*ndim,ndir
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);

      ! Set indices and directions
      idim1=mod(idir,ndir)+1
      idim2=mod(idir+1,ndir)+1

      jxC^L=ixC^L+kr(idim1,^D);
      ixCp^L=ixC^L+kr(idim2,^D);

      ! Interpolate sqrt gamma to the edges

      sqrtg(ixC^S)=1.d0

   !   select case(idim1)
   ! { case(^D)
   !     sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
   !       block%mSurface^D %sqrtgamma(ixC^S)+&
   !       block%mSurface^D %sqrtgamma(ixCp^S))
   ! \}
   !   end select
   !   
   !   select case(idim2)
   ! { case(^D)
   !     sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
   !       block%mSurface^D %sqrtgamma(ixC^S)+&
   !       block%mSurface^D %sqrtgamma(jxC^S))
   ! \}
   !   end select

      ! Reconstruct transverse transport velocities
      call reconstruct(ixI^L,ixC^L,idim2,vbarC(ixI^S,idim1,1),&
               vtilL(ixI^S,2),vtilR(ixI^S,2))

      call reconstruct(ixI^L,ixC^L,idim1,vbarC(ixI^S,idim2,2),&
               vtilL(ixI^S,1),vtilR(ixI^S,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      call reconstruct(ixI^L,ixC^L,idim2,bfaces(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))

      call reconstruct(ixI^L,ixC^L,idim1,bfaces(ixI^S,idim2),&
               btilL(ixI^S,idim2),btilR(ixI^S,idim2))

      ! Take the maximum characteristic

      cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
      cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

      cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
      cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
     
      ! Calculate eletric field
      if (idir <= ndim) then
         dxidir = dxlevel(idir)
      else
         dxidir = 1.0d0
      end if

      fE(ixC^S,idir)=qdt * dxidir * &
                   sqrtg(ixC^S) * (&
                   -(cp(ixC^S,1)*vtilL(ixC^S,1)*btilL(ixC^S,idim2) &
                   + cm(ixC^S,1)*vtilR(ixC^S,1)*btilR(ixC^S,idim2) &
                   -cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))&
                   /(cp(ixC^S,1) + cm(ixC^S,1)) &
                   +(cp(ixC^S,2)*vtilL(ixC^S,2)*btilL(ixC^S,idim1) &
                   + cm(ixC^S,2)*vtilR(ixC^S,2)*btilR(ixC^S,idim1) &
                   -cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))&
                   /(cp(ixC^S,2) + cm(ixC^S,2)) &
                   )

      if (.not.slab) then
        where(abs(x(ixC^S,r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixC^S,idir)=zero
        end where
      end if

    end do

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face: interal(fE dot dl)

    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=1,ndir ! Direction of line integral
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate
  end subroutine update_faces_uct2

  !> update faces UCT2 av
  subroutine updatefacesuct2av(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fC,fE,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qdt
    type(state)                        :: s
    double precision, intent(in)       :: vbarC(ixI^S,1:ndir,2)
    double precision, intent(in)       :: cbarmin(ixI^S,ndim)
    double precision, intent(in)       :: cbarmax(ixI^S,ndim)
    double precision, intent(in)       :: fC(ixI^S,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixI^S,1:ndir)

    double precision                   :: vtilL(ixI^S,2)
    double precision                   :: vtilR(ixI^S,2)
    double precision                   :: btilL(s%ixGs^S,ndim)
    double precision                   :: btilR(s%ixGs^S,ndim)
    double precision                   :: sqrtg(s%ixGs^S)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: fEc(ixI^S)
    double precision                   :: fEvar(ixI^S)
    double precision                   :: fEmin(ixI^S)
    double precision                   :: fEmax(ixI^S)
    integer                            :: ixGs^L
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,jxCp^L,ixCm^L
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixI^S,1:ndim), dxidir
    integer                            :: i^D

    associate(bfaces=>s%ws,w=>s%w,x=>s%x)
    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixGs^L=s%ixGs^L;

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    fE(ixI^S,1:ndir)=zero

    do idir=7-2*ndim,ndir
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);

      ! Set indices and directions
      idim1=mod(idir,ndir)+1
      idim2=mod(idir+1,ndir)+1

      jxC^L=ixC^L+kr(idim1,^D);
      ixCp^L=ixC^L+kr(idim2,^D);
      jxCp^L=jxC^L+kr(idim2,^D);

      ! Interpolate sqrt gamma to the edges

      sqrtg(ixC^S)=1.d0

  !    select case(idim1)
  !  { case(^D)
  !      sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
  !        block%mSurface^D %sqrtgamma(ixC^S)+&
  !        block%mSurface^D %sqrtgamma(ixCp^S))
  !  \}
  !    end select
  !    
  !    select case(idim2)
  !  { case(^D)
  !      sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
  !        block%mSurface^D %sqrtgamma(ixC^S)+&
  !        block%mSurface^D %sqrtgamma(jxC^S))
  !  \}
  !    end select

      ! Reconstruct transverse transport velocities
      call reconstruct(ixI^L,ixC^L,idim2,vbarC(ixI^S,idim1,1),&
               vtilL(ixI^S,2),vtilR(ixI^S,2))

      call reconstruct(ixI^L,ixC^L,idim1,vbarC(ixI^S,idim2,2),&
               vtilL(ixI^S,1),vtilR(ixI^S,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      call reconstruct(ixI^L,ixC^L,idim2,bfaces(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))

      call reconstruct(ixI^L,ixC^L,idim1,bfaces(ixI^S,idim2),&
               btilL(ixI^S,idim2),btilR(ixI^S,idim2))

      ! Take the maximum characteristic

      cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1))
      cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

      cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2))
      cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))
     
      ! Calculate elctric field
      if (idir <= ndim) then
         dxidir = dxlevel(idir)
      else
         dxidir = 1.0d0
      end if

      fE(ixC^S,idir)=qdt * dxidir * &
                   sqrtg(ixC^S) * (&
                   -(cp(ixC^S,1)*vtilL(ixC^S,1)*btilL(ixC^S,idim2) &
                   + cm(ixC^S,1)*vtilR(ixC^S,1)*btilR(ixC^S,idim2) &
                   -cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))&
                   /(cp(ixC^S,1) + cm(ixC^S,1)) &
                   +(cp(ixC^S,2)*vtilL(ixC^S,2)*btilL(ixC^S,idim1) &
                   + cm(ixC^S,2)*vtilR(ixC^S,2)*btilR(ixC^S,idim1) &
                   -cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))&
                   /(cp(ixC^S,2) + cm(ixC^S,2)) &
                   )

      
      if (.not.slab) then
        where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
          fE(ixC^S,idir)=zero
        end where
      end if

      ! Allowed range for the electric field in that direction.
      ! For simplicity, use twice cell-centred values, to allow
      ! some overshooting due to reconstruction
      ! Recall: 1,2,3 --> idir, idim1, idim2

   !   fEc(ixI^S) = myM%beta(idim2)%elem(ixI^S)*w(ixI^S,iw_mag(idim1)) &
   !              - myM%beta(idim1)%elem(ixI^S)*w(ixI^S,iw_mag(idim2))

   !   fEvar(ixI^S) = 2.0*max(abs(w(ixI^S,iw_mag(idim2)) - w(ixI^S,iw_mag(idim1))),&
   !                      abs(w(ixI^S,iw_mag(idim1))),&
   !                      abs(w(ixI^S,iw_mag(idim2))),1e-6)

   !   fEmax(ixI^S) = fEc(ixI^S) + myM%alpha(ixI^S)*fEvar(ixI^S)

   !   fEmin(ixI^S) = fEc(ixI^S) - myM%alpha(ixI^S)*fEvar(ixI^S)

   !   ! For the edge-allocated electric fields, take the maximum and
   !   ! minimum of the surrounding cells

   ! {#IFNDEF D1

   !   fEmax(ixC^S) = qdt * dxidir * sqrtg(ixC^S) * max(fEmax(ixC^S),fEmax(jxC^S),fEmax(ixCp^S),fEmax(jxCp^S)) 

   !   fEmin(ixC^S) = qdt * dxidir * sqrtg(ixC^S) * min(fEmin(ixC^S),fEmin(jxC^S),fEmin(ixCp^S),fEmin(jxCp^S)) 

   ! }
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Fall back to Balsara-Spicer averaging
      ! in case of problems
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (any(.not.(fE(ixC^S,idir).lt.fEmax(ixC^S).and.fE(ixC^S,idir).gt.fEmin(ixC^S)))) then
    !     write(*,*) mype, 'WARNING: Had to fall back to Balsara-Spicer','it=',it
         iwdim1 = iw_mag(idim1); iwdim2 = iw_mag(idim2)
         where (.not.(fE(ixC^S,idir).lt.fEmax(ixC^S).and.fE(ixC^S,idir).gt.fEmin(ixC^S)))
            fE(ixC^S,idir)=qdt*quarter*((fC(ixC^S,iwdim1,idim2)&
                 +fC(jxC^S,iwdim1,idim2))/dxlevel(idim1)&
                 -(fC(ixC^S,iwdim2,idim1)+fC(ixCp^S,iwdim2,idim1))&
                 /dxlevel(idim2))
         end where
      end if

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Before clipping the electric field,
      ! forcefully remove NaNs that could be there
      ! Brutalistic (probably never helps
      ! as NaNs will also be in other variables then)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (any(fE(ixC^S,idir) .ne. fE(ixC^S,idir))) then
    !     write(*,*) mype, 'WARNING: Had to forecefully set E=0','it=',it
         where (fE(ixC^S,idir) .ne. fE(ixC^S,idir))
            fE(ixC^S,idir) = zero
         end where
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If even Balsara-Spicer yields electric
      ! field out of bounds, reset to a fraction
      ! of the allowed field
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (any(.not.(fE(ixC^S,idir).lt.fEmax(ixC^S).and.fE(ixC^S,idir).gt.fEmin(ixC^S)))) then
    !     write(*,*) mype, 'WARNING: Clipping the electric field','it=',it
         where (fE(ixC^S,idir).gt.fEmax(ixC^S))
            fE(ixC^S,idir)=fEmin(ixC^S) + 0.8d0 * (fEmax(ixC^S) - fEmin(ixC^S))
         end where

         where (fE(ixC^S,idir).lt.fEmin(ixC^S))
            fE(ixC^S,idir)=fEmin(ixC^S) + 0.2d0 * (fEmax(ixC^S) - fEmin(ixC^S))
         end where

      end if

    end do

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=1,ndir ! Direction of line integral
           
          ! Assemble indices
          hxC^L=ixC^L-kr(idim2,^D);
          ! Add line integrals in direction idir
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)&
                           *(fE(ixC^S,idir)&
                            -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
       ixCmax^D=ixOmax^D;
       ixCmin^D=ixOmin^D-kr(idim1,^D);
       where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
         circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
       elsewhere
         circ(ixC^S,idim1)=zero
       end where
       ! Time update
        bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

  end associate
  end subroutine updatefacesuct2av

  !> update faces UCT1
  subroutine updatefacesuct1(ixI^L,ixO^L,qdt,vbarRC,vbarLC,cbarmin,cbarmax,fE,s)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: qdt
    type(state)                        :: s
    double precision, intent(in)       :: vbarRC(ixI^S,1:ndir,2)
    double precision, intent(in)       :: vbarLC(ixI^S,1:ndir,2)
    double precision, intent(in)       :: cbarmin(ixI^S,ndim)
    double precision, intent(in)       :: cbarmax(ixI^S,ndim)
    double precision, intent(inout)    :: fE(ixI^S,1:ndir)

    double precision                   :: vtilRL(ixI^S,2)
    double precision                   :: vtilLL(ixI^S,2)
    double precision                   :: vtilRR(ixI^S,2)
    double precision                   :: vtilLR(ixI^S,2)
    double precision                   :: btilL(s%ixGs^S,ndim)
    double precision                   :: btilR(s%ixGs^S,ndim)
    double precision                   :: sqrtg(s%ixGs^S)
    double precision                   :: cp(ixI^S,2)
    double precision                   :: cm(ixI^S,2)
    double precision                   :: ELL(ixI^S)
    double precision                   :: ELR(ixI^S)
    double precision                   :: ERL(ixI^S)
    double precision                   :: ERR(ixI^S)
    integer                            :: ixGs^L
    integer                            :: hxC^L,ixC^L,ixCp^L,jxC^L,ixCm^L
    integer                            :: idim1,idim2,idir
    double precision                   :: circ(ixI^S,1:ndim), dxidir

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ixGs^L=s%ixGs^L;

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    fE(ixI^S,1:ndir)=zero

    do idir=7-2*ndim,ndir
      ! Indices
      ! idir: electric field component
      ! idim1: where the first reconstruction was done
      ! idim2: where we will do the second reconstruction
      ! cyclic permutation: idim1,idim2,idir=1,2,3

      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-1+kr(idir,^D);

     ! Set indices and directions
      idim1=mod(idir,ndir)+1
      idim2=mod(idir+1,ndir)+1

      jxC^L=ixC^L+kr(idim1,^D);
      ixCp^L=ixC^L+kr(idim2,^D);

      ! Interpolate sqrt gamma to the edges

      sqrtg(ixC^S)=1.d0

   !   select case(idim1)
   ! { case(^D)
   !     sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
   !       block%mSurface^D %sqrtgamma(ixC^S)+&
   !       block%mSurface^D %sqrtgamma(ixCp^S))
   ! \}
   !   end select
      
   !   select case(idim2)
   ! { case(^D)
   !     sqrtg(ixC^S)=sqrtg(ixC^S)+quarter*(&
   !       block%mSurface^D %sqrtgamma(ixC^S)+&
   !       block%mSurface^D %sqrtgamma(jxC^S))
   ! \}
   !   end select

      ! Reconstruct velocities
      call reconstruct(ixI^L,ixC^L,idim2,vbarLC(ixI^S,idir,1),&
               vtilLL(ixI^S,1),vtilLR(ixI^S,1)) 
      call reconstruct(ixI^L,ixC^L,idim2,vbarRC(ixI^S,idir,1),&
               vtilRL(ixI^S,1),vtilRR(ixI^S,1))

      call reconstruct(ixI^L,ixC^L,idim2,vbarLC(ixI^S,idir,2),&
               vtilLL(ixI^S,2),vtilLR(ixI^S,2)) 
      call reconstruct(ixI^L,ixC^L,idim2,vbarRC(ixI^S,idir,2),&
               vtilRL(ixI^S,2),vtilRR(ixI^S,2))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      call reconstruct(ixI^L,ixC^L,idim2,bfaces(ixI^S,idim1),&
               btilL(ixI^S,idim1),btilR(ixI^S,idim1))
      call reconstruct(ixI^L,ixC^L,idim1,bfaces(ixI^S,idim2),&
               btilL(ixI^S,idim2),btilR(ixI^S,idim2))

      ! Take the maximum characteristic
      cm(ixC^S,1)=max(cbarmin(ixCp^S,idim1),cbarmin(ixC^S,idim1)) 
      cp(ixC^S,1)=max(cbarmax(ixCp^S,idim1),cbarmax(ixC^S,idim1))

      cm(ixC^S,2)=max(cbarmin(jxC^S,idim2),cbarmin(ixC^S,idim2)) 
      cp(ixC^S,2)=max(cbarmax(jxC^S,idim2),cbarmax(ixC^S,idim2))

      ! Calculate component idir of vxB partial
      ELL(ixC^S)= vtilLL(ixC^S,2)*btilL(ixC^S,idim1)&
           -vtilLL(ixC^S,1)*btilL(ixC^S,idim2)

      ELR(ixC^S)= vtilLR(ixC^S,2)*btilR(ixC^S,idim1)&
           -vtilLR(ixC^S,1)*btilL(ixC^S,idim2)

      ERL(ixC^S)= vtilRL(ixC^S,2)*btilL(ixC^S,idim1)&
           -vtilRL(ixC^S,1)*btilR(ixC^S,idim2)

      ERR(ixC^S)= vtilRR(ixC^S,2)*btilR(ixC^S,idim1)&
           -vtilRR(ixC^S,1)*btilR(ixC^S,idim2)

      ! For 3D, use interval
      if (idir <= ndim) then
         dxidir = dxlevel(idir)
      else
         dxidir = 1.0d0
      end if
      
      ! Calculate elctric field
      fE(ixC^S,idir)=qdt * dxidir * &
       sqrtg(ixC^S) * (&
       (cp(ixC^S,1)*cp(ixC^S,2)*ELL(ixC^S)&
       +cp(ixC^S,1)*cm(ixC^S,2)*ELR(ixC^S)&
       +cm(ixC^S,1)*cp(ixC^S,2)*ERL(ixC^S)&
       +cm(ixC^S,1)*cm(ixC^S,2)*ERR(ixC^S)&
       )/&
        ((cp(ixC^S,1)+cm(ixC^S,1))*&
         (cp(ixC^S,2)+cm(ixC^S,2)))&
       +(cp(ixC^S,1)*cm(ixC^S,1)*(btilR(ixC^S,idim2)-btilL(ixC^S,idim2)))/&
         (cp(ixC^S,1) + cm(ixC^S,1))&
       -(cp(ixC^S,2)*cm(ixC^S,2)*(btilR(ixC^S,idim1)-btilL(ixC^S,idim1)))/&
         (cp(ixC^S,2) + cm(ixC^S,2)))

      if (.not.slab) then
        where(abs(x(ixC^S,r_)+half*dxlevel(r_))<1.0d-9)
          fE(ixC^S,idir)=zero
        end where
      end if

    end do

    circ(ixI^S,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face 
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      do idim2=1,ndim
        do idir=1,ndir ! Direction of line integral
         ! Assemble indices
         hxC^L=ixC^L-kr(idim2,^D);
         ! Add line integrals in direction idir
         circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                          +lvc(idim1,idim2,idir)&
                          *(fE(ixC^S,idir)&
                           -fE(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(s%surfaceC(ixC^S,idim1) > 1.0d-9*s%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/s%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ! Time update
      bfaces(ixC^S,idim1)=bfaces(ixC^S,idim1)-circ(ixC^S,idim1)
    end do

    end associate
  end subroutine updatefacesuct1

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
    double precision                   :: circ(ixI^S,1:ndim), dxidir(ixI^S)

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
          if (idir <= ndim) then
            dxidir(ixC^S) = block%dsC(ixC^S,idir)
          else
            dxidir(ixC^S) = 1.0d0
          end if
          circ(ixC^S,idim1)=circ(ixC^S,idim1)&
                           +lvc(idim1,idim2,idir)*dxidir(ixC^S) &
                           *(A(ixC^S,idir)&
                            -A(hxC^S,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get B
    do idim1=1,ndim
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(idim1,^D);
      where(block%surfaceC(ixC^S,idim1) > 1.0d-9*block%dvolume(ixC^S))
        circ(ixC^S,idim1)=circ(ixC^S,idim1)/block%surfaceC(ixC^S,idim1)
      elsewhere
        circ(ixC^S,idim1)=zero
      end where
      ws(ixC^S,idim1) = circ(ixC^S,idim1)
    end do

  end subroutine b_from_vectorpotentialA

  !> Routine for easily printing an array, just for debugging, erase later
  subroutine printarray(ixI^L,ixO^L,array)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: array(ixI^S)
    integer                       :: i,j,k

    {^IFTWOD
    do j=ixOmax2,ixOmin2,-1
    }
       do i=ixOmin1,ixOmax1
    {^IFONED
          write (*,'(I2,E17.9,$)') i,array(i)
    }
    {^IFTWOD
          write (*,'(I2,I2,E17.9,$)') i,j,array(i,j)
    }
    {^IFTHREED
          do k=ixOmin3,ixOmax3
            write (*,'(I2,I2,I2,E17.9,$)') i,j,k,array(i,j,k)
          end do
          print *, ' '
          }
       end do

       print *, '---'
    {^IFTWOD
    end do
    }

    write (*,'(A)') '------------------------------------------------------------'
  end subroutine printarray

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
