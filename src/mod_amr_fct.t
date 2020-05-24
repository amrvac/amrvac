module mod_amr_fct
  implicit none
  private

  type facealloc
    double precision, dimension(:^D&), pointer :: face
  end type facealloc

  type fake_neighbors
    integer :: igrid
    integer :: ipe
  end type fake_neighbors

  type(facealloc), dimension(:,:,:), allocatable, public :: pface

  type(fake_neighbors), dimension(:^D&,:,:), allocatable, public :: fine_neighbors

  integer, dimension(:,:^D&,:), allocatable, public :: old_neighbor

  integer :: itag, isend, irecv
  integer :: nrecv, nsend, ibuf_recv, ibuf_send, ibuf_send_next
  integer, dimension(^ND) :: isize
  integer, dimension(:), allocatable :: recvrequest, sendrequest
  integer, dimension(:,:), allocatable :: recvstatus, sendstatus
  double precision, allocatable :: recvbuffer(:), sendbuffer(:)

  public :: store_faces
  public :: comm_faces
  public :: end_comm_faces
  public :: deallocateBfaces
  public :: old_neighbors
  public :: prolong_2nd_stg
  public :: already_fine

contains
  !> This subroutine performs a 2nd order prolongation for a staggered field F,
  !> preserving the divergence of the coarse cell.
  !> This is useful for preserving DivF=0.
  !> If DivF=f(x), a different algorithm must be used.
  subroutine prolong_2nd_stg(sCo,sFi,ixCo^Lin,ixFi^Lin,dxCo^D,xComin^D,dxFi^D,xFimin^D,ghost,fine_^Lin)
    use mod_global_parameters
    use mod_physics

    logical, intent(in)          :: ghost
    integer, intent(in)          :: ixCo^Lin, ixFi^Lin
    double precision, intent(in) :: dxCo^D, xComin^D, dxFi^D, xFimin^D
    type(state), intent(in)      :: sCo
    type(state), intent(inout)   :: sFi

    logical, optional :: fine_^Lin
    logical           :: fine_^L

    double precision :: eta^D, invdxCo^D
    integer :: ixCo^L,ixFi^L
    integer :: idim1,idim2,ix^DE,idim3,ixFis^L,ixGs^L,ixCos^L,ixFisC^L
    integer :: ixCosV^L(1:ndim),ixFisV^L(1:ndim)
    integer :: hxCos^L,jxCos^L,ixCosE^L,ixFisE^L,hxFisC^L,jxFisC^L,ipxFisC^L,ixCosC^L,imxFisC^L,jpxFisC^L,jmxFisC^L,hpxFisC^L
    integer :: hxFi^L,jxFi^L,hijxFi^L,hjixFi^L,hjjxFi^L
    integer :: iihxFi^L,iijxFi^L,ijhxFi^L,ijjxFi^L,ihixFi^L,ijixFi^L,ihjxFi^L
    integer :: jihxFi^L,jijxFi^L,jjhxFi^L,jjjxFi^L,jhixFi^L,jjixFi^L,jhjxFi^L
    double precision :: bfluxCo(sCo%ixGs^S,nws),bfluxFi(sFi%ixGs^S,nws)
    double precision :: slopes(sCo%ixGs^S,ndim),B_energy_change(ixG^T)
    {^IFTHREED
    ! Directional bias, see pdf
    ! These arrays should have ranges ixO^S
    ! Since they are still not known, we give ranges ixG^T,
    ! even though it is wasteful
    double precision :: sigmau(ixG^T),sigmad(ixG^T),sigma(ixG^T,1:ndim), alpha(ixG^T,1:ndim)
    ! Auxiliary arrays for magnetic fluxes
    double precision :: F1(ixG^T),F2(ixG^T),F3(ixG^T),F4(ixG^T)
    }

    {^IFONED
    call mpistop("CT prolongation not implemented in 1D. But CT is not needed.")
    }

    {^NOONED
    ! Note on the indices:
    ! ixCo  Cells where 
    !       divergence-preserving prolongation will be applied.
    ! ixFi  Fine cells which correspond to that extent.
    ! ixCoE For 'expanded', faces that need to be used to calculate
    !       the slopes for the prolongation
    ! ixCos
    ! ixFis
    ! ixCosV And
    ! ixFisV For 'variable', since the ranges depend on the direction of the faces
    ! ixCosC For 'copy',
    ! ixFisC For 'copy', to fill the information from the fine grid needed in the internal faces prolongation step.
    ! At the end, ixFisC is used to copy the assign the magnetic
    ! fields components to the state structure.

    fine_^L=.false.;

    {if(present(fine_min^Din)) fine_min^D=fine_min^Din;}
    {if(present(fine_max^Din)) fine_max^D=fine_max^Din;}

    ! When filling ghost cells, ixFi^L are given.
    ! When refining the block, ixCo^L are given.
    if (ghost) then
      ixFi^L=ixFi^Lin;
      {ixComin^D=int((ixFimin^D+nghostcells+1)/2);}
      {ixComax^D=int((ixFimax^D+nghostcells+1)/2);}
    else
      ixCo^L=ixCo^Lin;
      ixFi^L=ixM^LL;
    end if

    ! Expanded range for staggered variables
    ixCosmin^D=ixComin^D-1;
    ixCosmax^D=ixComax^D;

    ixFismin^D=ixFimin^D-1;
    ixFismax^D=ixFimax^D;


    associate(wCos=>sCo%ws, wFis=>sFi%ws,wCo=>sCo%w, wFi=>sFi%w)
    ! Assemble general indices
    ixGsmin^D=sFi%ixGsmin^D;
    ixGsmax^D=sFi%ixGsmax^D;

    do idim1=1,ndim
      ixCosVmin^D(idim1)=ixComin^D-kr(^D,idim1);
      ixCosVmax^D(idim1)=ixComax^D;
      ixFisVmin^D(idim1)=ixFimin^D-kr(^D,idim1);
      ixFisVmax^D(idim1)=ixFimax^D;
    end do

    ! Initialize auxiliary arrays at zero
    bfluxCo = zero
    bfluxFi = zero 
    slopes  = zero 


    invdxCo^D=1.d0/dxCo^D;

    ! Fill coarse magnetic flux array
    do idim1=1,ndim
      ! Fill information from parts already at the fine level
      ! First set up indices
      if (ghost) then
        ixFisEmin^D=max(1-kr(^D,idim1),ixFisVmin^D(idim1)-2*(1-kr(^D,idim1)));
        ixFisEmax^D=min(ixGsmax^D,ixFisVmax^D(idim1)+2*(1-kr(^D,idim1)));
        ixCosE^L=int((ixFisE^L+nghostcells+1)/2);
      else
        ixCosE^L=ixCosV^L(idim1)^LADD(1-kr(idim1,^D));
        ixFisE^L=ixFisV^L(idim1)^LADD2*(1-kr(idim1,^D));
      end if
      ! Convert fine fields to fluxes
      bfluxFi(ixFisE^S,idim1)=wFis(ixFisE^S,idim1)*sFi%surfaceC(ixFisE^S,idim1)
      {^IFTWOD
      idim2=1+mod(idim1,2)
      }
      {^IFTHREED
      idim2=1+mod(idim1,3)
      idim3=1+mod(idim1+1,3)
      }
      bfluxCo(ixCosE^S,idim1) = zero
      ! Add fine fluxes sharing the same fine face
     {do ix^DE=0,1\}
        ixFisC^L=ixFisE^L+ix2*kr(idim2,^D);
        {^IFTHREED 
        ixFisC^L=ixFisC^L+ix3*kr(idim3,^D);
        }
        bfluxCo(ixCosE^S,idim1)=bfluxCo(ixCosE^S,idim1)+bfluxFi(ixFisCmin^D:ixFisCmax^D:2,idim1)
     {end do^DE&\}
    end do

    ! Omit indices for already refined face, if any
    {
    if (fine_min^D) then
      ixCosVmin^D(^D)=ixCosVmin^D(^D)+1
      ixFisVmin^D(^D)=ixFisVmin^D(^D)+2
    end if
    if (fine_max^D) then
      ixCosVmax^D(^D)=ixCosVmax^D(^D)-1
      ixFisVmax^D(^D)=ixFisVmax^D(^D)-2
    end if
    \}

    do idim1=1,ndim
      ixCosE^L=ixCosV^L(idim1);
      ! Omit part already refined
      {
      if (^D/=idim1) then
        if ((.not.fine_min^D).or.(.not.ghost)) then
         ixCosEmin^D=ixCosVmin^D(idim1)-1
        end if
        if ((.not.fine_max^D).or.(.not.ghost)) then
         ixCosEmax^D=ixCosVmax^D(idim1)+1
        end if
      end if
      \}
      ! Fill coarse flux array from coarse field
      bfluxCo(ixCosE^S,idim1)=wCos(ixCosE^S,idim1)*sCo%surfaceC(ixCosE^S,idim1)
    end do
    ! Finished filling coarse flux array

    ! Distribute coarse fluxes among fine fluxes
    ! There are too many loops here, perhaps optimize later
    do idim1=1,ndim
       ixCos^L=ixCosV^L(idim1);
       do idim2=1,ndim
         if(idim1==idim2) cycle
         ! Calculate slope in direction idim2
         ! Set up indices
         jxCos^L=ixCos^L+kr(idim2,^D);
         hxCos^L=ixCos^L-kr(idim2,^D);
         slopes(ixCos^S,idim2)=0.125d0*(bfluxCo(jxCos^S,idim1)-bfluxCo(hxCos^S,idim1))
    {^IFTWOD
         do ix2=0,1
           ixFisCmin^D=ixFisVmin^D(idim1)+ix2*kr(^D,idim2);
           ixFisCmax^D=ixFisVmax^D(idim1)+ix2*kr(^D,idim2);
           bfluxFi(ixFisCmin^D:ixFisCmax^D:2,idim1)=half*(bfluxCo(ixCos^S,idim1)&
              +(2*ix2-1)*slopes(ixCos^S,idim2))
         end do
    }
    {^IFTHREED
         do idim3=1,ndim
           ! Calculate slope in direction idim3
           ! Set up indices
           jxCos^L=ixCos^L+kr(idim3,^D);
           hxCos^L=ixCos^L-kr(idim3,^D);
           slopes(ixCos^S,idim3)=0.125d0*(bfluxCo(jxCos^S,idim1)-bfluxCo(hxCos^S,idim1))
           if(lvc(idim1,idim2,idim3)<1) cycle
           do ix2=0,1
             do ix3=0,1
              {ixFisCmin^D=ixFisVmin^D(idim1)+ix2*kr(^D,idim2)+ix3*kr(^D,idim3);}
              {ixFisCmax^D=ixFisVmax^D(idim1)+ix2*kr(^D,idim2)+ix3*kr(^D,idim3);}
               bfluxFi(ixFisCmin^D:ixFisCmax^D:2,idim1)=quarter*bfluxCo(ixCos^S,idim1)&
                 +quarter*(2*ix2-1)*slopes(ixCos^S,idim2)&
                 +quarter*(2*ix3-1)*slopes(ixCos^S,idim3)
             end do
           end do
         end do
    }
       end do
    end do

    ! Calculate interior fine fluxes
    {^IFTHREED
    ! Directional bias for nonlinear prolongation
    do idim1=1,ndim
      do idim2=1,ndim
        do idim3=1,ndim
          if (lvc(idim1,idim2,idim3)<1) cycle
          ! Set up indices
          hxFi^L=ixFi^L-kr(idim1,^D);
          jxFi^L=ixFi^L+kr(idim1,^D);

          hijxFi^L=hxFi^L+kr(idim3,^D);
          hjixFi^L=hxFi^L+kr(idim2,^D);
          hjjxFi^L=hijxFi^L+kr(idim2,^D);

          iihxFi^L=ixFi^L-kr(idim3,^D);
          iijxFi^L=ixFi^L+kr(idim3,^D);
          ihixFi^L=ixFi^L-kr(idim2,^D);
          ihjxFi^L=ihixFi^L+kr(idim3,^D);
          ijixFi^L=ixFi^L+kr(idim2,^D);
          ijhxFi^L=ijixFi^L-kr(idim3,^D);
          ijjxFi^L=ijixFi^L+kr(idim3,^D);

          jihxFi^L=jxFi^L-kr(idim3,^D);
          jijxFi^L=jxFi^L+kr(idim3,^D);
          jhixFi^L=jxFi^L-kr(idim2,^D);
          jhjxFi^L=jhixFi^L+kr(idim3,^D);
          jjixFi^L=jxFi^L+kr(idim2,^D);
          jjhxFi^L=jjixFi^L-kr(idim3,^D);
          jjjxFi^L=jjixFi^L+kr(idim3,^D);

          sigmau(ixCo^S)=&
              abs(bfluxFi(jhixFimin^D:jhixFimax^D:2,idim2))&
            + abs(bfluxFi(jhjxFimin^D:jhjxFimax^D:2,idim2))&
            + abs(bfluxFi(jjixFimin^D:jjixFimax^D:2,idim2))&
            + abs(bfluxFi(jjjxFimin^D:jjjxFimax^D:2,idim2))&
            + abs(bfluxFi(jihxFimin^D:jihxFimax^D:2,idim3))&
            + abs(bfluxFi(jijxFimin^D:jijxFimax^D:2,idim3))&
            + abs(bfluxFi(jjhxFimin^D:jjhxFimax^D:2,idim3))&
            + abs(bfluxFi(jjjxFimin^D:jjjxFimax^D:2,idim3))

          sigmad(ixCo^S)=&
              abs(bfluxFi(ihixFimin^D:ihixFimax^D:2,idim2))&
            + abs(bfluxFi(ihjxFimin^D:ihjxFimax^D:2,idim2))&
            + abs(bfluxFi(ijixFimin^D:ijixFimax^D:2,idim2))&
            + abs(bfluxFi(ijjxFimin^D:ijjxFimax^D:2,idim2))&
            + abs(bfluxFi(iihxFimin^D:iihxFimax^D:2,idim3))&
            + abs(bfluxFi(iijxFimin^D:iijxFimax^D:2,idim3))&
            + abs(bfluxFi(ijhxFimin^D:ijhxFimax^D:2,idim3))&
            + abs(bfluxFi(ijjxFimin^D:ijjxFimax^D:2,idim3))

          sigma(ixCo^S,idim1)=sigmau(ixCo^S)+sigmad(ixCo^S)
          where(sigma(ixCo^S,idim1)/=zero)
            sigma(ixCo^S,idim1)=abs(sigmau(ixCo^S)-sigmad(ixCo^S))/sigma(ixCo^S,idim1)
          elsewhere
            sigma(ixCo^S,idim1)=zero
          end where

        end do
      end do
    end do

    ! Directional bias for Toth-Roe prolongation
    do idim1=1,ndim
      do idim2=1,ndim
        do idim3=1,ndim
          if (lvc(idim1,idim2,idim3)<1) cycle
    !       Nonlinear
    !        alpha(ixCo^S,idim1)=sigma(ixCo^S,idim2)-sigma(ixCo^S,idim3)
    !       Homogeneous
    !        alpha(ixCo^S,idim1)=0.d0
    !       Toth-Roe
            alpha(ixCo^S,idim1)=(sCo%dx(ixCo^S,idim2)-sCo%dx(ixCo^S,idim3))/(sCo%dx(ixCo^S,idim2)+sCo%dx(ixCo^S,idim3))
        end do
      end do
    end do
    }

    do idim1=1,ndim
      do idim2=1,ndim
        {^IFTWOD
        if (idim1==idim2) cycle
        ixFisCmin^D=ixFismin^D+1;
        ixFisCmax^D=ixFismax^D-1;
        jxFisC^L=ixFisC^L+kr(idim1,^D);
        hxFisC^L=ixFisC^L-kr(idim1,^D);
        ipxFisC^L=ixFisC^L+kr(idim2,^D);
        imxFisC^L=ixFisC^L-kr(idim2,^D);
        jpxFisC^L=jxFisC^L+kr(idim2,^D);
        jmxFisC^L=jxFisC^L-kr(idim2,^D);
        hpxFisC^L=hxFisC^L+kr(idim2,^D);

        bfluxFi(ixFisCmin^D:ixFisCmax^D:2,idim1)=&
         half*(bfluxFi(jxFisCmin^D:jxFisCmax^D:2,idim1)&
              +bfluxFi(hxFisCmin^D:hxFisCmax^D:2,idim1))&
         -quarter*(bfluxFi(ipxFisCmin^D:ipxFisCmax^D:2,idim2)&
                  -bfluxFi(jpxFisCmin^D:jpxFisCmax^D:2,idim2)&
                  -bfluxFi(imxFisCmin^D:imxFisCmax^D:2,idim2)&
                  +bfluxFi(jmxFisCmin^D:jmxFisCmax^D:2,idim2))

        bfluxFi(ipxFisCmin^D:ipxFisCmax^D:2,idim1)=&
         half*(bfluxFi(jpxFisCmin^D:jpxFisCmax^D:2,idim1)&
              +bfluxFi(hpxFisCmin^D:hpxFisCmax^D:2,idim1))&
         -quarter*(bfluxFi(ipxFisCmin^D:ipxFisCmax^D:2,idim2)&
                  -bfluxFi(jpxFisCmin^D:jpxFisCmax^D:2,idim2)&
                  -bfluxFi(imxFisCmin^D:imxFisCmax^D:2,idim2)&
                  +bfluxFi(jmxFisCmin^D:jmxFisCmax^D:2,idim2))
        }
        {^IFTHREED
        do idim3=1,ndim
          if (lvc(idim1,idim2,idim3)<1) cycle
          ! Set up indices
          hxFi^L=ixFi^L-kr(idim1,^D);
          jxFi^L=ixFi^L+kr(idim1,^D);

          hijxFi^L=hxFi^L+kr(idim3,^D);
          hjixFi^L=hxFi^L+kr(idim2,^D);
          hjjxFi^L=hijxFi^L+kr(idim2,^D);

          iihxFi^L=ixFi^L-kr(idim3,^D);
          iijxFi^L=ixFi^L+kr(idim3,^D);
          ihixFi^L=ixFi^L-kr(idim2,^D);
          ihjxFi^L=ihixFi^L+kr(idim3,^D);
          ijixFi^L=ixFi^L+kr(idim2,^D);
          ijhxFi^L=ijixFi^L-kr(idim3,^D);
          ijjxFi^L=ijixFi^L+kr(idim3,^D);

          jihxFi^L=jxFi^L-kr(idim3,^D);
          jijxFi^L=jxFi^L+kr(idim3,^D);
          jhixFi^L=jxFi^L-kr(idim2,^D);
          jhjxFi^L=jhixFi^L+kr(idim3,^D);
          jjixFi^L=jxFi^L+kr(idim2,^D);
          jjhxFi^L=jjixFi^L-kr(idim3,^D);
          jjjxFi^L=jjixFi^L+kr(idim3,^D);

          ! Prolongation formulas
          F1(ixCo^S)=bfluxFi(ihixFimin^D:ihixFimax^D:2,idim2)&
                    -bfluxFi(jhixFimin^D:jhixFimax^D:2,idim2)&
                    -bfluxFi(ijixFimin^D:ijixFimax^D:2,idim2)&
                    +bfluxFi(jjixFimin^D:jjixFimax^D:2,idim2)

          F2(ixCo^S)=bfluxFi(ihjxFimin^D:ihjxFimax^D:2,idim2)&
                    -bfluxFi(jhjxFimin^D:jhjxFimax^D:2,idim2)&
                    -bfluxFi(ijjxFimin^D:ijjxFimax^D:2,idim2)&
                    +bfluxFi(jjjxFimin^D:jjjxFimax^D:2,idim2)

          F3(ixCo^S)=bfluxFi(iihxFimin^D:iihxFimax^D:2,idim3)&
                    -bfluxFi(jihxFimin^D:jihxFimax^D:2,idim3)&
                    -bfluxFi(iijxFimin^D:iijxFimax^D:2,idim3)&
                    +bfluxFi(jijxFimin^D:jijxFimax^D:2,idim3)

          F4(ixCo^S)=bfluxFi(ijhxFimin^D:ijhxFimax^D:2,idim3)&
                    -bfluxFi(jjhxFimin^D:jjhxFimax^D:2,idim3)&
                    -bfluxFi(ijjxFimin^D:ijjxFimax^D:2,idim3)&
                    +bfluxFi(jjjxFimin^D:jjjxFimax^D:2,idim3)

          bfluxFi(ixFimin^D:ixFimax^D:2,idim1)=&
             half*(bfluxFi(hxFimin^D:hxFimax^D:2,idim1)+bfluxFi(jxFimin^D:jxFimax^D:2,idim1))&
             +6.25d-2*((3.d0+alpha(ixCo^S,idim2))*F1(ixCo^S)&
                      +(1.d0-alpha(ixCo^S,idim2))*F2(ixCo^S)&
                      +(3.d0-alpha(ixCo^S,idim3))*F3(ixCo^S)&
                      +(1.d0+alpha(ixCo^S,idim3))*F4(ixCo^S))

          bfluxFi(ijixFimin^D:ijixFimax^D:2,idim1)=&
             half*(bfluxFi(hjixFimin^D:hjixFimax^D:2,idim1)+bfluxFi(jjixFimin^D:jjixFimax^D:2,idim1))&
             +6.25d-2*((3.d0+alpha(ixCo^S,idim2))*F1(ixCo^S)&
                      +(1.d0-alpha(ixCo^S,idim2))*F2(ixCo^S)&
                      +(1.d0+alpha(ixCo^S,idim3))*F3(ixCo^S)&
                      +(3.d0-alpha(ixCo^S,idim3))*F4(ixCo^S))

          bfluxFi(iijxFimin^D:iijxFimax^D:2,idim1)=&
             half*(bfluxFi(hijxFimin^D:hijxFimax^D:2,idim1)+bfluxFi(jijxFimin^D:jijxFimax^D:2,idim1))&
             +6.25d-2*((1.d0-alpha(ixCo^S,idim2))*F1(ixCo^S)&
                      +(3.d0+alpha(ixCo^S,idim2))*F2(ixCo^S)&
                      +(3.d0-alpha(ixCo^S,idim3))*F3(ixCo^S)&
                      +(1.d0+alpha(ixCo^S,idim3))*F4(ixCo^S))

          bfluxFi(ijjxFimin^D:ijjxFimax^D:2,idim1)=&
             half*(bfluxFi(hjjxFimin^D:hjjxFimax^D:2,idim1)+bfluxFi(jjjxFimin^D:jjjxFimax^D:2,idim1))&
             +6.25d-2*((1.d0-alpha(ixCo^S,idim2))*F1(ixCo^S)&
                      +(3.d0+alpha(ixCo^S,idim2))*F2(ixCo^S)&
                      +(1.d0+alpha(ixCo^S,idim3))*F3(ixCo^S)&
                      +(3.d0-alpha(ixCo^S,idim3))*F4(ixCo^S))
        end do
        }
      end do
    end do

    ! Go back to magnetic fields
    do idim1=1,ndim
      ixFisCmax^D=ixFimax^D;
      ixFisCmin^D=ixFimin^D-kr(^D,idim1);
      where(sFi%surfaceC(ixFisC^S,idim1)/=zero)
        wFis(ixFisC^S,idim1)=bfluxFi(ixFisC^S,idim1)/sFi%surfaceC(ixFisC^S,idim1)
      elsewhere
        wFis(ixFisC^S,idim1)=zero
      end where
    end do

    if(phys_energy.and. .not.prolongprimitive) then
      B_energy_change(ixFi^S)=0.5d0*sum(wFi(ixFi^S,iw_mag(:))**2,dim=ndim+1)
    end if
    call phys_face_to_center(ixFi^L,sFi)
    if(phys_energy.and. .not.prolongprimitive) then
      B_energy_change(ixFi^S)=0.5d0*sum(wFi(ixFi^S,iw_mag(:))**2,dim=ndim+1)-&
        B_energy_change(ixFi^S)
      wFi(ixFi^S,iw_e)=wFi(ixFi^S,iw_e)+B_energy_change(ixFi^S)
    end if

    end associate

    } ! END NOONED
  end subroutine prolong_2nd_stg

  !> To achive consistency and thus conservation of divergence,
  !> when refining a block we take into account the faces of the
  !> already fine neighbours, if any. This routine stores them.
  subroutine store_faces
    use mod_forest, only: refine 
    use mod_global_parameters
    integer :: igrid, iigrid, idims, iside, ineighbor, ipe_neighbor
    integer :: nx^D, i^D, ic^D, inc^D

    if (npe>1) then
      nsend_fc=0
      nrecv_fc=0
    end if

    ! Size of the block face
    nx^D=ixMhi^D-ixMlo^D+1;

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ! Check whether it is necessary to store any block face, i.e.
      ! if any coarser neighbour is going to be refined 
     {do iside=1,2
        i^DD=kr(^DD,^D)*(2*iside-3);
        if (neighbor_pole(i^DD,igrid)/=0) cycle
        if (neighbor_type(i^DD,igrid)==neighbor_coarse) then
          ineighbor   =neighbor(1,i^DD,igrid)
          ipe_neighbor=neighbor(2,i^DD,igrid)
          if (refine(ineighbor,ipe_neighbor)) then
            allocate(pface(iside,^D,igrid)%face(1^D%1:nx^DD))
            !! Store the faces
            if (iside==1) then !! left side
              pface(iside,^D,igrid)%face(1^D%1:nx^DD)=&
                ps(igrid)%ws(ixMlo^D-1^D%ixM^T,^D)
            else !! right side
              pface(iside,^D,igrid)%face(1^D%1:nx^DD)=&
                ps(igrid)%ws(ixMhi^D^D%ixM^T,^D)
            end if
            if (ipe_neighbor/=mype) nsend_fc(^D)=nsend_fc(^D)+1
          end if
        end if
      end do\}

      ! If a grid is going to be refined,
      ! remember what are its neighbours.
      if (refine(igrid,mype)) then
        ! Initialize neighbour array
        fine_neighbors(:^D&,:,igrid)%igrid=-1
        fine_neighbors(:^D&,:,igrid)%ipe=-1
        do idims=1,ndim
          do iside=1,2
            i^D=kr(^D,idims)*(2*iside-3);
            if (neighbor_pole(i^D,igrid)/=0) cycle
            if (neighbor_type(i^D,igrid)==neighbor_fine) then
             {do ic^DB=1+int((1+i^DB)/2),2-int((1-i^DB)/2)
                inc^DB=ic^DB+i^DB\}
                ineighbor=neighbor_child(1,inc^D,igrid)
                ipe_neighbor=neighbor_child(2,inc^D,igrid)

                fine_neighbors(ic^D,idims,igrid)%igrid= ineighbor
                fine_neighbors(ic^D,idims,igrid)%ipe=ipe_neighbor

                if (ipe_neighbor/=mype) nrecv_fc(idims)=nrecv_fc(idims)+1
             {end do\}
            end if
          end do
        end do
      end if

    end do

  end subroutine store_faces

  !> When refining a coarse block with fine neighbours, it is necessary
  !> prolong consistently with the already fine faces.
  !> This routine takes care of the communication of such faces.
  subroutine comm_faces
    use mod_forest, only: refine
    use mod_global_parameters

    integer                   :: iigrid,igrid,ineighbor,ipe_neighbor
    integer                   :: idims,iside,i^D,ic^D,inc^D,nx^D
    integer                   :: recvsize, sendsize

    ! Communicate the block faces to achieve consistency when refining
    ! Initialize communication structures

    nrecv=0
    nsend=0
    recvsize=0
    sendsize=0

    do idims=1,ndim
       select case (idims)
       {case (^D)
          nrecv=nrecv+nrecv_fc(^D)
          nsend=nsend+nsend_fc(^D)
          nx^D=1^D%nx^DD=ixMhi^DD-ixMlo^DD+1;
          isize(^D)={nx^DD*}
          recvsize=recvsize+nrecv_fc(^D)*isize(^D)
          sendsize=sendsize+nsend_fc(^D)*isize(^D)
       \}
       end select
    end do

    if (nrecv>0) then
    ! Allocate receive buffer
      allocate(recvbuffer(recvsize),recvstatus(MPI_STATUS_SIZE,nrecv), &
               recvrequest(nrecv))
      recvrequest=MPI_REQUEST_NULL
      ibuf_recv=1
      irecv=0

    ! Receive
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        if (refine(igrid,mype)) then
         {do ic^DB=1,2\}
            ! Only one of the sides will be necessary,
            ! so we do the loop only over dimensions, instead of
            ! over dimensions and sizes as in the routines
            ! old_neighbors and already_fine.
            do idims=1,ndim
              ipe_neighbor=fine_neighbors(ic^D,idims,igrid)%ipe
              ineighbor   =fine_neighbors(ic^D,idims,igrid)%igrid
              if (ineighbor>0.and.ipe_neighbor/=mype) then
               {if (idims==^D) iside=ic^D\}
                !!! Check indices
                i^D=kr(^D,idims)*(2*iside-3);
                if (neighbor_pole(i^D,igrid)/=0) cycle
                inc^D=ic^D+i^D;
                irecv=irecv+1
                itag=4**^ND*(igrid-1)+{inc^D*4**(^D-1)+}

                call MPI_IRECV(recvbuffer(ibuf_recv),isize(idims), &
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                               icomm,recvrequest(irecv),ierrmpi)
                ibuf_recv=ibuf_recv+isize(idims)
              end if
            end do
         {end do\}
        end if
      end do

    end if

    if (nsend>0) then
    ! Allocate send buffer
      allocate(sendbuffer(sendsize),sendstatus(MPI_STATUS_SIZE,nsend),sendrequest(nsend))
      sendrequest=MPI_REQUEST_NULL
      isend=0
      ibuf_send=1
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ! Check whether it is necessary to store any block face, i.e.
        ! if any coarser neighbour is going to be refined 
       {do iside=1,2
          i^DD=kr(^DD,^D)*(2*iside-3);
          ! When there is a pole, faces are always zero and this is not necessary
          if (neighbor_pole(i^DD,igrid)/=0) cycle
          if (neighbor_type(i^DD,igrid)==neighbor_coarse) then
            ineighbor   =neighbor(1,i^DD,igrid)
            ipe_neighbor=neighbor(2,i^DD,igrid)
            if (refine(ineighbor,ipe_neighbor)) then
              if (ipe_neighbor/=mype) then
                ic^DD=1+modulo(node(pig^DD_,igrid)-1,2);
                inc^D=-2*i^D+ic^D^D%inc^DD=ic^DD;
                itag=4**^ND*(ineighbor-1)+{inc^DD*4**(^DD-1)+}
                isend=isend+1
                ibuf_send_next=ibuf_send+isize(^D)
                sendbuffer(ibuf_send:ibuf_send_next-1)=&
                reshape(pface(iside,^D,igrid)%face,(/isize(^D)/))
                call MPI_ISEND(sendbuffer(ibuf_send),isize(^D), &
                               MPI_DOUBLE_PRECISION,ipe_neighbor,itag, &
                               icomm,sendrequest(isend),ierrmpi)
                ibuf_send=ibuf_send_next
              end if
            end if
          end if
        end do\}
      end do
    end if

    ! Waitalls
    if (nrecv>0) then
       call MPI_WAITALL(nrecv,recvrequest,recvstatus,ierrmpi)
       deallocate(recvstatus,recvrequest)
       ibuf_recv=1
    end if

    if (nsend>0) then
       call MPI_WAITALL(nsend,sendrequest,sendstatus,ierrmpi)
       deallocate(sendbuffer,sendstatus,sendrequest)
    end if

  end subroutine comm_faces

  subroutine end_comm_faces
    use mod_global_parameters
    ! Deallocate receive buffer
    if (nrecv>0) deallocate(recvbuffer)
  end subroutine end_comm_faces

  subroutine deallocateBfaces
    use mod_global_parameters
    integer :: igrid, iigrid, iside

    do iigrid=1,igridstail; igrid=igrids(iigrid);
     {do iside=1,2
        if (associated(pface(iside,^D,igrid)%face)) then
          deallocate(pface(iside,^D,igrid)%face)
        end if
      end do\}
    end do

  end subroutine deallocateBfaces

  subroutine old_neighbors(child_igrid,child_ipe,igrid,ipe)
    use mod_global_parameters
    integer, dimension(2^D&), intent(in) :: child_igrid, child_ipe
    integer, intent(in) :: igrid, ipe
    integer :: iside, i^D, ic^D

    {do ic^DB=1,2\}
      old_neighbor(:,:^D&,child_igrid(ic^D))=-1
     {do iside=1,2
        if (ic^D==iside) then
          i^DD=kr(^DD,^D)*(2*iside-3);
          old_neighbor(1,i^DD,child_igrid(ic^DD))=&
             fine_neighbors(ic^DD,^D,igrid)%igrid
          old_neighbor(2,i^DD,child_igrid(ic^DD))=&
             fine_neighbors(ic^DD,^D,igrid)%ipe
        end if
      end do\}
    {end do\}

  end subroutine old_neighbors

  !> This routine fills the fine faces before prolonging.
  !> It is the face equivalent of fix_conserve 
  subroutine already_fine(sFi,ichild,fine_^L)
    use mod_forest
    use mod_global_parameters
    type(tree_node_ptr) :: tree
    type(state) :: sFi
    integer, intent(in) :: ichild
    logical :: fine_^L

    integer :: ineighbor,ipe_neighbor,ibufnext
    integer :: iside,iotherside,i^D,nx^D

    ! Size of the block face
    nx^D=ixMhi^D-ixMlo^D+1;

    ! Initialise everything to zero and false
    fine_min^D=.false.;
    fine_max^D=.false.;
    sFi%ws=zero

    {^NOONED
    ! This face communication is not needed in 1D
   {do iside=1,2
      i^DD=kr(^DD,^D)*(2*iside-3);
      ! This is not necessary at the pole.
      ! We are inside a loop over the children, so the grid index is ichild
      if (neighbor_pole(i^DD,ichild)/=0) cycle

      ! Get old ipe and igrid of neighbour from array fake_neighbor
      ! Then plug it into the structure pfaces and get the faces 
      ineighbor   =old_neighbor(1,i^DD,ichild)
      ipe_neighbor=old_neighbor(2,i^DD,ichild)

      iotherside=3-iside
      if (ineighbor>0) then
        if (iside==1) then ! left side
          if (ipe_neighbor==mype) then
            sFi%ws(ixMlo^D-1^D%ixM^T,^D)=pface(iotherside,^D,ineighbor)%face(1^D%1:nx^DD)
          else
            ibufnext=ibuf_recv+isize(^D)
            sFi%ws(ixMlo^D-1^D%ixM^T,^D)=reshape(&
                 source=recvbuffer(ibuf_recv:ibufnext-1),&
                 shape=shape(sFi%ws(ixMlo^D-1^D%ixM^T,^D)))
            ibuf_recv=ibufnext
          end if
          fine_min^D=.true.
        else ! right side
          if (ipe_neighbor==mype) then
            sFi%ws(ixMhi^D^D%ixM^T,^D)=pface(iotherside,^D,ineighbor)%face(1^D%1:nx^DD)
          else
            ibufnext=ibuf_recv+isize(^D)
            sFi%ws(ixMhi^D^D%ixM^T,^D)=reshape(&
                 source=recvbuffer(ibuf_recv:ibufnext-1),&
                 shape=shape(sFi%ws(ixMhi^D^D%ixM^T,^D)))
            ibuf_recv=ibufnext
          end if
          fine_max^D=.true.
        end if
      end if
    end do\}
    }
  end subroutine already_fine

end module mod_amr_fct
