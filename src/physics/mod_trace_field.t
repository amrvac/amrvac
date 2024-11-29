module mod_trace_field
  use mod_global_parameters
  use mod_physics
  implicit none

contains

  subroutine trace_field_multi(xfm,wPm,wLm,dL,numL,numP,nwP,nwL,forwardm,ftype,tcondi)
    ! trace multiple field lines
    ! xfm: locations of points at the field lines. User should provide xfm(1:numL,1,1:ndim)
    !   as seed points, then field line wills be traced from the seed points. xfm(1:numL,1,:) 
    !   are user given and other points are given by the subroutine as feedback.
    ! numL: number of field lines user wants to trace. user given
    ! numP: maximum number of points at the field line. User defined. Note that not every
    !   point of numP is valid, as the tracing will stop when the field line leave the
    !   simulation box. The number of valid point is given in wLm(1)
    ! wPm: point variables, which have values at each point of xfm. The way to get wP is
    !   user defined, with help of IO given by subroutine set_field_w in mod_usr_methods.
    !   User can calculate the density/temperature at the point and then store the values
    !   into wPm, or do something else.
    ! nwP: number of point variables. user given
    ! wLm: line variables, variables for the lines rather than for each point of the lines.
    !   For example, wLm(1) stores the number of valid points at the field lines. The way to 
    !   get wLm is also user defined with set_field_w, the same as wPm. User can calculate
    !   the maximum density/temperature at the field lines and stores it in wLm
    ! nwL: number of line variables. user given
    ! dL: length step for field line tracing. user given
    ! forwardm: true--trace field forward; false--trace field line backward. user given
    ! ftype: type of field user wants to trace. User can trace velocity field by setting
    !   ftype='Vfield' or trace magnetic field by setting ftype='Bfield'. It is possible 
    !   to trace other fields, e.g. electric field, where user can define the field with
    !   IO given by subroutine set_field in mod_usr_methods. user given
    ! tcondi: user given
    use mod_particle_base

    integer, intent(in) :: numL,numP,nwP,nwL
    double precision, intent(inout) :: xfm(numL,numP,ndim),wPm(numL,numP,nwP),wLm(numL,1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forwardm(numL)
    character(len=std_len), intent(in) :: ftype,tcondi

    integer :: indomain,ipe_now,igrid_now,igrid,j,iL
    integer :: ipoint_in,ipoint_out,numSend,nRT,nRTmax
    double precision :: x3d(3),statusF(4+ndim),statusL(numL,4+ndim+nwL),statusS(numL,4+ndim+nwL)
    logical :: continueL(numL),myL(numL)
    logical :: stopT,forward
    integer :: ipointm(numL),igridm(numL)
    double precision :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, allocatable :: data_send(:,:,:),data_recv(:,:,:)

    if (tcondi/='TRAC') then
      wPm=zero
    else
      wPm=-1
    endif
    wLm=zero
    xfm(1:numL,2:numP,:)=zero
    stopT=.TRUE.
    myL=.FALSE.
    xf=zero
    wP=zero

    ! find the pe and igrid for the first point
    do iL=1,numL
      indomain=0
      wLm(iL,1)=0
      {if (xfm(iL,1,^DB)>=xprobmin^DB .and. xfm(iL,1,^DB)<xprobmax^DB) indomain=indomain+1\}
      if (indomain==ndim) then
        if (tcondi/='TRAC') wLm(iL,1)=1
        continueL(iL)=.TRUE.
        ! find pe and igrid
        x3d=0.d0
        do j=1,ndim
          x3d(j)=xfm(iL,1,j)
        enddo
        call find_particle_ipe(x3d,igrid_now,ipe_now)
        stopT=.FALSE.
        ipointm(iL)=1
        igridm(iL)=igrid_now
        if (mype==ipe_now) then 
          myL(iL)=.TRUE.
        else
          xfm(iL,1,:)=zero
        endif
      else
        continueL(iL)=.FALSE.
        wLm(iL,1)=zero
      endif
    enddo

    do while(stopT .eqv. .FALSE.)
      ! tracing multiple field lines inside pe
      statusS=zero
      do iL=1,numL
        if (myL(iL) .and. continueL(iL)) then
          igrid=igridm(iL)
          ipoint_in=ipointm(iL)
          xf(ipoint_in,:)=xfm(iL,ipoint_in,:)
          wL(:)=wLm(iL,:)
          forward=forwardm(iL)
          statusF=zero
          call find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)
          ipoint_out=int(statusF(1))
          xfm(iL,ipoint_in:ipoint_out-1,:)=xf(ipoint_in:ipoint_out-1,:)
          wPm(iL,ipoint_in:ipoint_out-1,:)=wP(ipoint_in:ipoint_out-1,:)
          ! status for each field line
          ! 1: index of next point
          ! 2: ipe of next point
          ! 3: igrid of next point
          ! 4: 1-> continue tracing; 0-> stop tracing
          ! 5:4+ndim: coordinate of next point
          ! 4+ndim+1:4+ndim+nwL: wL(2:1+nwL)
          ! for TRAC nwL=2 -> wL(2): current Tcoff; wL(3): Tmax
          ! for TRAC nwP=2 -> wP(:,1):ipe; wP(:,2):igrid
          statusS(iL,1:4+ndim)=statusF(1:4+ndim)
          statusS(iL,4+ndim+1:4+ndim+nwL)=wL(2:1+nwL)
          if (tcondi=='TRAC') wLm(iL,1)=ipoint_out-1
        endif
      enddo

      ! comunicating tracing results
      numSend=numL*(4+ndim+nwL)
      call MPI_ALLREDUCE(statusS,statusL,numSend,MPI_DOUBLE_PRECISION,&
                       MPI_SUM,icomm,ierrmpi)

      ! for next step
      stopT=.TRUE.
      myL=.FALSE.
      do iL=1,numL
        if (continueL(iL)) then
          ipointm(iL)=int(statusL(iL,1))
          if (mype==int(statusL(iL,2))) myL(iL)=.TRUE.
          igridm(iL)=int(statusL(iL,3))
          if (int(statusL(iL,4))==0) then
            stopT=.FALSE.
          else
            continueL(iL)=.FALSE.
          endif
          if (myL(iL)) xfm(iL,ipointm(iL),1:ndim)=statusL(iL,4+1:4+ndim)
          if (tcondi/='TRAC') wLm(iL,1)=ipointm(iL)-1
          wLm(iL,2:1+nwL)=statusL(iL,4+ndim+1:4+ndim+nwL)
        endif
      enddo
    enddo

    ! communication after tracing
    if (tcondi/='TRAC') then
      nRTmax=0
      do iL=1,numL
        if (nRTmax<int(wLm(iL,1))) nRTmax=int(wLm(iL,1))
      enddo
      numSend=numL*nRTmax*(ndim+nwP)

      allocate(data_send(numL,nRTmax,ndim+nwP),data_recv(numL,nRTmax,ndim+nwP))
      data_send(:,:,:)=zero
      do iL=1,numL
        nRT=int(wLm(iL,1))
        data_send(iL,1:nRT,1:ndim)=xfm(iL,1:nRT,1:ndim)
        if (nwP>0) data_send(iL,1:nRT,1+ndim:ndim+nwP)=wPm(iL,1:nRT,1:nwP)
      enddo
      call MPI_ALLREDUCE(data_send,data_recv,numSend,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      do iL=1,numL
        nRT=int(wLm(iL,1))
        xfm(iL,1:nRT,1:ndim)=data_recv(iL,1:nRT,1:ndim)
        if (nwP>0) wPm(iL,1:nRT,1:nwP)=data_recv(iL,1:nRT,1+ndim:ndim+nwP)
      enddo
      deallocate(data_send,data_recv)
    endif

  end subroutine trace_field_multi

  subroutine trace_field_single(xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi)
    ! trace a field line
    ! xf: locations of points at the field line. User should provide xf(1,1:ndim)
    !   as seed point, then field line will be traced from the seed point. xf(1,:) is 
    !   user given and xf(2:wL(1),:) are given by the subroutine as feedback.
    ! numP: maximum number of points at the field line. User defined. Note that not every
    !   point of numP is valid, as the tracing will stop when the field line leave the
    !   simulation box. The number of valid point is given in wL(1)
    ! wP: point variables, which have values at each point of xf. The way to get wP is
    !   user defined, with help of IO given by subroutine set_field_w in mod_usr_methods.
    !   User can calculate the density/temperature at the point and then store the values
    !   into wP, or do something else.
    ! nwP: number of point variables. user given
    ! wL: line variables, variables for the line rather than for each point of the line.
    !   For example, wL(1) stores the number of valid points at the field line. The way to 
    !   get wL is also user defined with set_field_w, the same as wP. User can calculate
    !   the maximum density/temperature at the field line and stores it in wL
    ! nwL: number of line variables. user given
    ! dL: length step for field line tracing. user given
    ! forward: true--trace field forward; false--trace field line backward. user given
    ! ftype: type of field user wants to trace. User can trace velocity field by setting
    !   ftype='Vfield' or trace magnetic field by setting ftype='Bfield'. It is possible 
    !   to trace other fields, e.g. electric field, where user can define the field with
    !   IO given by subroutine set_field in mod_usr_methods. user given
    ! tcondi: user given
    use mod_usr_methods
    use mod_particle_base

    integer, intent(in) :: numP,nwP,nwL
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi

    integer :: indomain,ipoint_in,ipe_now,igrid_now,igrid,j
    integer :: ipoint_out,ipe_next,igrid_next,numRT
    double precision :: x3d(3),statusF(4+ndim),status_bcast(4+ndim+nwL)
    logical :: stopT
    double precision, allocatable :: data_send(:,:),data_recv(:,:)

    wP=zero
    wL=zero
    xf(2:numP,:)=zero

    ! check whether or the first point is inside simulation box. if yes, find
    ! the pe and igrid for the point
    indomain=0
    wL(1)=0
    {if (xf(1,^DB)>=xprobmin^DB .and. xf(1,^DB)<xprobmax^DB) indomain=indomain+1\}
    if (indomain==ndim) then
      wL(1)=1

       ! find pe and igrid
       x3d=0.d0
       do j=1,ndim
         x3d(j)=xf(1,j)
       enddo
      call find_particle_ipe(x3d,igrid_now,ipe_now)
      stopT=.FALSE.
      ipoint_in=1
      if (mype/=ipe_now) xf(1,:)=zero
    else
      if (mype==0) then
        call MPISTOP('Field tracing error: given point is not in simulation box!')
      endif
    endif


    ! other points in field line
    do while(stopT .eqv. .FALSE.)

      if (mype==ipe_now) then
        igrid=igrid_now
        ! looking for points in one pe
        call find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)
        status_bcast(1:4+ndim)=statusF(1:4+ndim)
        status_bcast(4+ndim+1:4+ndim+nwL)=wL(2:1+nwL)
      endif
      ! comunication
      call MPI_BCAST(status_bcast,4+ndim+nwL,MPI_DOUBLE_PRECISION,ipe_now,icomm,ierrmpi)
      statusF(1:4+ndim)=status_bcast(1:4+ndim)
      wL(2:1+nwL)=status_bcast(4+ndim+1:4+ndim+nwL)

      ! prepare for next step
      ipoint_out=int(statusF(1))
      ipe_next=int(statusF(2))
      igrid_next=int(statusF(3))
      if (int(statusF(4))==1) then
        stopT=.TRUE.
        wL(1)=ipoint_out-1
      endif
      if (mype==ipe_next) then
        do j=1,ndim
          xf(ipoint_out,j)=statusF(4+j)
        enddo
      else
        xf(ipoint_out,:)=zero
      endif

      ! pe and grid of next point
      ipe_now=ipe_next
      igrid_now=igrid_next
      ipoint_in=ipoint_out
    enddo

    if (tcondi/='TRAC') then
      numRT=int(wL(1))
      allocate(data_send(numRT,ndim+nwP),data_recv(numRT,ndim+nwP))
      data_send(:,:)=zero
      data_recv(:,:)=zero
      data_send(1:numRT,1:ndim)=xf(1:numRT,1:ndim)
      if (nwP>0) data_send(1:numRT,1+ndim:ndim+nwP)=wP(1:numRT,1:nwP)
      call MPI_ALLREDUCE(data_send,data_recv,numRT*(ndim+nwP),MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
      xf(1:numRT,1:ndim)=data_recv(1:numRT,1:ndim)
      if (nwP>0) wP(1:numRT,1:nwP)=data_recv(1:numRT,1+ndim:ndim+nwP)
      deallocate(data_send,data_recv)
    endif

  end subroutine trace_field_single

  subroutine find_points_in_pe(igrid,ipoint_in,xf,wP,wL,dL,numP,nwP,nwL,forward,ftype,tcondi,statusF)

    integer, intent(inout) :: igrid
    integer, intent(in) :: ipoint_in,numP,nwP,nwL
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi
    double precision, intent(inout) :: statusF(4+ndim)

    integer :: ipe_next,igrid_next,ip_in,ip_out,j,indomain
    logical :: newpe,stopT
    double precision :: xfout(ndim)

    ip_in=ipoint_in
    newpe=.FALSE.

    do while(newpe .eqv. .FALSE.)
      ! looking for points in given grid    
      call find_points_interp(igrid,ip_in,ip_out,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
      ip_in=ip_out

      ! when next point is out of given grid, find next grid  
      indomain=0
      {if (xf(ip_out,^DB)>=xprobmin^DB .and. xf(ip_out,^DB)<=xprobmax^DB) indomain=indomain+1\}      
      if (ip_out<numP .and. indomain==ndim) then
        if (tcondi/='TRAC') then
          stopT=.FALSE.
          xfout=xf(ip_out,:)
          call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopT)
        else
          if (xf(ip_out,ndim)>phys_trac_mask) then
            newpe=.TRUE.
            stopT=.TRUE.
          else
            stopT=.FALSE.
            xfout=xf(ip_out,:)
            call find_next_grid(igrid,igrid_next,ipe_next,xfout,newpe,stopT)
          endif
        endif
      else
        newpe=.TRUE.
        stopT=.TRUE.
      endif

      if (newpe) then
        statusF(1)=ip_out
        statusF(2)=ipe_next
        statusF(3)=igrid_next
        statusF(4)=0
        if (stopT) statusF(4)=1
        do j=1,ndim
          statusF(4+j)=xf(ip_out,j)
        enddo
      endif

      if (newpe .eqv. .FALSE.) igrid=igrid_next
    enddo

  end subroutine find_points_in_pe

  subroutine find_next_grid(igrid,igrid_next,ipe_next,xf1,newpe,stopT)
    ! check the grid and pe of next point
    use mod_usr_methods
    use mod_global_parameters
    use mod_forest

    integer, intent(inout) :: igrid,igrid_next,ipe_next
    double precision, intent(in) :: xf1(ndim)
    logical, intent(inout) :: newpe,stopT

    double precision :: dxb^D,xb^L,xbmid^D
    integer :: idn^D,my_neighbor_type,inblock
    integer :: ic^D,inc^D,ipe_neighbor,igrid_neighbor
    double precision :: xbn^L

    ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
    ^D&xbmax^D=rnode(rpxmax^D_,igrid)\
    inblock=0

    ! direction of next grid
    idn^D=0\
    {if (xf1(^D)<=xbmin^D) idn^D=-1\}
    {if (xf1(^D)>=xbmax^D) idn^D=1\}
    my_neighbor_type=neighbor_type(idn^D,igrid)
    igrid_neighbor=neighbor(1,idn^D,igrid)
    ipe_neighbor=neighbor(2,idn^D,igrid)

    ! ipe and igrid of next grid
    select case(my_neighbor_type)
    case (neighbor_boundary)
      ! next point is not in simulation box
      newpe=.TRUE.
      stopT=.TRUE.

    case(neighbor_coarse)
      ! neighbor grid has lower refinement level      
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_sibling)
      ! neighbor grid has lower refinement level 
      igrid_next=igrid_neighbor
      ipe_next=ipe_neighbor
      if (mype==ipe_neighbor) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif

    case(neighbor_fine)
      ! neighbor grid has higher refinement level 
      {xbmid^D=(xbmin^D+xbmax^D)/2.d0\}
      ^D&inc^D=1\
      {if (xf1(^D)<=xbmin^D) inc^D=0\}
      {if (xf1(^D)>xbmin^D .and. xf1(^D)<=xbmid^D) inc^D=1\}
      {if (xf1(^D)>xbmid^D .and. xf1(^D)<xbmax^D) inc^D=2\}
      {if (xf1(^D)>=xbmax^D) inc^D=3\}
      ipe_next=neighbor_child(2,inc^D,igrid)
      igrid_next=neighbor_child(1,inc^D,igrid)
      if (mype==ipe_next) then
        newpe=.FALSE.
      else
        newpe=.TRUE.
      endif
    end select

  end subroutine find_next_grid

  subroutine find_points_interp(igrid,ip_in,ip_out,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in) :: igrid,ip_in,numP,nwP,nwL
    integer, intent(inout) :: ip_out
    double precision, intent(inout) :: xf(numP,ndim),wP(numP,nwP),wL(1+nwL)
    double precision, intent(in) :: dL
    logical, intent(in) :: forward
    character(len=std_len), intent(in) :: ftype,tcondi

    double precision :: dxb^D,xb^L
    integer          :: ip,inblock,ixI^L,ixO^L,j
    double precision :: field(ixg^T,ndir)
    double precision :: xs1(ndim),xs2(ndim),K1(ndim),K2(ndim)
    double precision :: xfpre(ndim),xfnow(ndim),xfnext(ndim)
    double precision :: Tpre,Tnow,Tnext,dTds,Lt,Lr,ds,T_bott,trac_delta

    ixI^L=ixG^LL;
    ixO^L=ixM^LL;
    ^D&dxb^D=rnode(rpdx^D_,igrid);
    ^D&xbmin^D=rnode(rpxmin^D_,igrid);
    ^D&xbmax^D=rnode(rpxmax^D_,igrid);

    if (tcondi/='TRAC') then
      ds=dL
    else
      ds=dxb^ND
      T_bott=2.d4/unit_temperature
      trac_delta=0.25d0
    endif

    ! main loop
    MAINLOOP: do ip=ip_in,numP-1

      ! integrate magnetic field with Runge-Kutta method
      xs1(:)=xf(ip,:)
      call get_K(xs1,igrid,K1,ixI^L,dxb^D,ftype)
      if (forward) then
        xs2(:)=xf(ip,:)+ds*K1(:)
      else
        xs2(:)=xf(ip,:)-ds*K1(:)
      endif
      call get_K(xs2,igrid,K2,ixI^L,dxb^D,ftype)
      if (forward) then
        xf(ip+1,:)=xf(ip,:)+ds*(0.5*K1(:)+0.5*K2(:))
      else
        xf(ip+1,:)=xf(ip,:)-ds*(0.5*K1(:)+0.5*K2(:))
      endif
      ip_out=ip+1

      ! get local values for variable via interpolation
      if (tcondi/='TRAC') then
        if (associated(usr_set_field_w)) then 
          call usr_set_field_w(igrid,ip,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,tcondi)
        endif
      else  ! get TRAC Tcoff
        wP(ip,1)=mype
        wP(ip,2)=igrid
        if (ip==ip_in) then
          if (forward) then
            xfpre(:)=xf(ip,:)-ds*K1(:)
          else
            xfpre(:)=xf(ip,:)+ds*K1(:)
          endif
          xfnow(:)=xf(ip,:)
          xfnext(:)=xf(ip+1,:)
          call get_T_loc_TRAC(igrid,xfpre,Tpre,ixI^L,dxb^D)
          call get_T_loc_TRAC(igrid,xfnow,Tnow,ixI^L,dxb^D)
          call get_T_loc_TRAC(igrid,xfnext,Tnext,ixI^L,dxb^D)
        else
          xfpre=xf(ip-1,:)
          xfnow(:)=xf(ip,:)
          xfnext(:)=xf(ip+1,:)
          Tpre=Tnow
          Tnow=Tnext
          call get_T_loc_TRAC(igrid,xfnext,Tnext,ixI^L,dxb^D)
        endif
        dTds=abs(Tnext-Tpre)/(2*ds)
        if (ip==1) then
          Lt=0.d0
          wL(2)=T_bott   ! current Tcofl
          wL(3)=Tnow     ! current Tlmax
        else
          Lt=0.d0
          if (dTds>0.d0) then
            Lt=Tnow/dTds
            Lr=ds
            ! renew cutoff temperature
            if(Lr>trac_delta*Lt) then
              if (Tnow>wL(2)) wL(2)=Tnow
            endif
          endif
          if (Tnow>wL(3)) wL(3)=Tnow
        endif
      endif

      ! exit loop if next point is not in this block
      inblock=0
      {if (xf(ip+1,^DB)>=xbmin^DB .and. xf(ip+1,^DB)<xbmax^DB) inblock=inblock+1\}
      if (tcondi=='TRAC' .and. xf(ip+1,ndim)>phys_trac_mask) inblock=0
      if (inblock/=ndim) exit MAINLOOP

    enddo MAINLOOP

  end subroutine find_points_interp

  subroutine get_K(xfn,igrid,K,ixI^L,dxb^D,ftype)
    use mod_usr_methods

    integer :: ixI^L,igrid
    double precision :: dxb^D
    double precision :: xfn(ndim),K(ndim)
    character(len=std_len) :: ftype

    integer          :: ixb^D,ix^D,ixbl^D,ixO^L,j
    double precision :: xd^D
    double precision :: field(0:1^D&,ndir),Fx(ndim),factor(0:1^D&)
    double precision :: vector(ixI^S,1:ndir)
    double precision :: Ftotal

    ^D&ixbl^D=floor((xfn(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(xfn(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;
    ^D&ixOmin^D=ixbl^D;
    ^D&ixOmax^D=ixbl^D+1;
    vector=zero

    field=zero
    if(ftype=='Bfield') then
      if(B0field) then
        if(allocated(iw_mag)) then
          vector(ixO^S,1:ndir)=ps(igrid)%w(ixO^S,iw_mag(1:ndir))+ps(igrid)%B0(ixO^S,1:ndir,0)
        else
          vector(ixO^S,1:ndir)=ps(igrid)%B0(ixO^S,1:ndir,0)
        endif
      else
        vector(ixO^S,1:ndir)=ps(igrid)%w(ixO^S,iw_mag(1:ndir))
      endif
    else if (ftype=='Vfield') then
      call phys_get_v(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,vector)
    endif
    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
      field(ix^D,1:ndir)=vector(ixbl^D+ix^D,1:ndir)
    {enddo\}    

    Fx=0.d0
    {do ix^DB=0,1\}
      do j=1,ndim
        Fx(j)=Fx(j)+field(ix^D,j)*factor(ix^D)
      enddo
    {enddo\}

    if (ftype=='Bfield' .or. ftype=='Vfield') then
      Fx=0.d0
      {do ix^DB=0,1\}
        do j=1,ndim
          Fx(j)=Fx(j)+field(ix^D,j)*factor(ix^D)
        enddo
      {enddo\}
    else if (associated(usr_set_field)) then
      call usr_set_field(xfn,igrid,Fx,ftype)
    else
      call MPISTOP('Field tracing error: wrong field type!')
    endif

    Ftotal=zero
    do j=1,ndim
      Ftotal=Ftotal+(Fx(j))**2
    enddo
    Ftotal=dsqrt(Ftotal)

    if (Ftotal==0.d0) then
      K=0.d0
      K(1)=1.d0
    else
      K(1:ndim)=Fx(1:ndim)/Ftotal
    endif

  end subroutine get_K

  subroutine get_T_loc_TRAC(igrid,xloc,Tloc,ixI^L,dxb^D)
    ! grid T has been calculated and stored in wextra(ixI^S,Tcoff_)
    ! see mod_trac
    integer, intent(in) :: igrid,ixI^L
    double precision, intent(inout) :: xloc(ndim)
    double precision, intent(inout) :: Tloc
    double precision, intent(in) :: dxb^D

    integer          :: ixb^D,ix^D,ixbl^D,j,ixO^L
    double precision :: xd^D
    double precision :: factor(0:1^D&),Tnear(0:1^D&)

    ^D&ixbl^D=floor((xloc(^D)-ps(igrid)%x(ixImin^DD,^D))/dxb^D)+ixImin^D;
    ^D&xd^D=(xloc(^D)-ps(igrid)%x(ixbl^DD,^D))/dxb^D;
    ^D&ixOmin^D=ixbl^D;
    ^D&ixOmax^D=ixOmin^D+1;

    {do ix^D=0,1\}
      factor(ix^D)={abs(1-ix^D-xd^D)*}
      Tnear(ix^D)=ps(igrid)%wextra(ixbl^D+ix^D,iw_tcoff)
    {enddo\}

    Tloc=0.d0
    ! interpolation
    {do ix^DB=0,1\}
      Tloc=Tloc+Tnear(ix^D)*factor(ix^D)
    {enddo\}

  end subroutine get_T_loc_TRAC

end module mod_trace_field
