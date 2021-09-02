!> module mod_magnetofriction.t
!> Purpose: use magnetofrictional method to relax 3D magnetic field to
!>          force-free field
!> 01.04.2016 developed by Chun Xia and Yang Guo
!> 04.10.2017 modulized by Chun Xia
!> Usage:
!> in amrvac.par:
!>   &methodlist
!>    time_stepper='onestep' ! time marching scheme, or 'twostep','threestep'
!>    flux_method=13*'cd4' ! or 'tvdlf', 'fd'
!>    limiter= 13*'koren' ! or 'vanleer','cada3','mp5' so on
!>   /
!>   &meshlist
!>    ditregrid=20 ! set iteration interval for adjusting AMR
!>   /
!>   &mhd_list
!>    mhd_magnetofriction=.true.
!>   /
!>   &mf_list
!>    mf_it_max=60000 ! set the maximum iteration number
!>    mf_ditsave=20000 ! set iteration interval for data output
!>    mf_cc=0.3    ! stability coefficient controls numerical stability
!>    mf_cy=0.2    ! frictional velocity coefficient
!>    mf_cdivb=0.1 ! divb cleaning coefficient controls diffusion speed of divb
!>   /
module mod_magnetofriction
  implicit none

  !> stability coefficient controls numerical stability
  double precision :: mf_cc
  !> frictional velocity coefficient
  double precision :: mf_cy, mf_cy_max
  !> divb cleaning coefficient controls diffusion speed of divb
  double precision :: mf_cdivb
  !> TVDLF dissipation coefficient controls the dissipation term
  double precision :: mf_tvdlfeps, mf_tvdlfeps_min
  !> time in magnetofriction process
  double precision :: tmf
  !> maximal speed for fd scheme
  double precision :: cmax_mype
  !> maximal speed for fd scheme
  double precision :: cmax_global

  !> Index of the density (in the w array)
  integer, private, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Indices of the magnetic field
  integer, allocatable, private, protected :: mag(:)

  integer :: mf_ditsave
  integer :: mf_it_max
  logical :: mf_advance
  logical :: fix_conserve_at_step = .true.

contains
  !> Read this module"s parameters from a file
  subroutine mf_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mf_list/ mf_ditsave, mf_it_max, mf_cc, mf_cy, mf_cy_max, mf_cdivb, mf_tvdlfeps_min

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mf_list, end=111)
111    close(unitpar)
    end do

  end subroutine mf_params_read

  !> Initialize the module
  subroutine magnetofriction_init()
    use mod_global_parameters
    use mod_usr_methods, only: usr_before_main_loop

    rho_=iw_rho
    allocate(mom(ndir))
    mom=iw_mom
    allocate(mag(ndir))
    mag=iw_mag

    mf_it_max=60000  ! set the maximum iteration number
    mf_ditsave=20000 ! set iteration interval for data output
    mf_cc=0.3d0      ! stability coefficient controls numerical stability
    mf_cy=0.2d0      ! frictional velocity coefficient. The default value is mf_cy itself.
    mf_cy_max=mf_cy  ! maximum of the frictional velocity coefficient
    mf_cdivb=0.1d0   ! divb cleaning coefficient controls diffusion speed of divb
    mf_tvdlfeps=1.d0 ! coefficient to control the TVDLF dissipation
    mf_tvdlfeps_min = mf_tvdlfeps ! minimum of the TVDLF dissipation coefficient

    call mf_params_read(par_files)

    usr_before_main_loop => magnetofriction

  end subroutine magnetofriction_init

  subroutine magnetofriction
    use mod_global_parameters
    use mod_physics
    use mod_ghostcells_update
    use mod_input_output

    double precision :: dvolume(ixG^T),dsurface(ixG^T),dvone
    double precision :: dtfff,dtfff_pe,dtnew,dx^D
    double precision :: cwsin_theta_new,cwsin_theta_old
    double precision :: sum_jbb,sum_jbb_ipe,sum_j,sum_j_ipe,sum_l_ipe,sum_l
    double precision :: f_i_ipe,f_i,volumepe,volume,tmpt,time_in
    double precision, external :: integral_grid
    integer :: i,iigrid, igrid, idims,ix^D,hxM^LL,fhmf,tmpit,i^D
    logical :: patchwi(ixG^T), stagger_flag

    ! not do fix conserve and getbc for staggered values if stagger is used
    stagger_flag=stagger_grid
    stagger_grid=.false.

    time_in=MPI_WTIME()
    if(mype==0) write(*,*) 'Evolving to force-free field using magnetofricitonal method...'
    if(prolongprimitive) call mpistop('use prolongprimitive=.false. in MF module')
    mf_advance=.false.
    dtfff=1.d-2
    tmpt=global_time
    tmpit=it
    tmf=global_time
    i=it
    if(snapshotini==0 .and. i==0) then
      call saveamrfile(1)
      call saveamrfile(2)
    end if
    mf_advance=.true.
    ! point bc mpi datatype to partial type for magnetic field
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    ! create bc mpi datatype for ghostcells update
    call create_bc_mpi_datatype(mag(1),ndir)
    ! point bc mpi datatype to partial type for velocity field
    type_send_srl=>type_send_srl_p2
    type_recv_srl=>type_recv_srl_p2
    type_send_r=>type_send_r_p2
    type_recv_r=>type_recv_r_p2
    type_send_p=>type_send_p_p2
    type_recv_p=>type_recv_p_p2
    ! create bc mpi datatype for ghostcells update
    call create_bc_mpi_datatype(mom(1),ndir)
    ! convert conservative variables to primitive ones which are used during MF
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call phys_to_primitive(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
    end do
    ! calculate magnetofrictional velocity
    call mf_velocity_update(dtfff)
    ! update velocity in ghost cells
    call getbc(tmf,0.d0,ps,mom(1),ndir)
    ! calculate initial values of metrics
    if(i==0) then
      call metrics
      call printlog_mf
    end if
    ! magnetofrictional loops
    do
      ! calculate time step based on Cmax= Alfven speed + abs(frictional speed)
      dtfff_pe=bigdouble
      cmax_mype=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        pso(igrid)%w(ixG^T,mag(:))=ps(igrid)%w(ixG^T,mag(:))
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        call getdtfff_courant(ps(igrid)%w,ps(igrid)%x,ixG^LL,ixM^LL,dtnew)
        dtfff_pe=min(dtfff_pe,dtnew)
      end do
      call MPI_ALLREDUCE(dtfff_pe,dtfff,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                         icomm,ierrmpi)
      call MPI_ALLREDUCE(cmax_mype,cmax_global,1,MPI_DOUBLE_PRECISION,&
                 MPI_MAX,icomm,ierrmpi)

      ! =======
      ! evolve
      ! =======
      call advectmf(1,ndim,tmf,dtfff)

      if(i>=10000)  then
        mf_tvdlfeps = 0.9998d0 * mf_tvdlfeps
        mf_tvdlfeps = max(mf_tvdlfeps_min,mf_tvdlfeps)
      end if
      if(i<=60000) then
        mf_cy=1.0001d0*mf_cy
        mf_cy = min(mf_cy_max,mf_cy)
      end if

      i=i+1
      tmf=tmf+dtfff
      if(mod(i,10)==0) then
        ! calculate metrics
        call metrics
        call printlog_mf
      end if
      if(mod(i,mf_ditsave)==0) then
        it=i
        global_time=tmf
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call phys_to_conserved(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do
        mf_advance=.false.
        call saveamrfile(1)
        call saveamrfile(2)
        do iigrid=1,igridstail; igrid=igrids(iigrid);
           call phys_to_primitive(ixG^LL,ixG^LL,ps(igrid)%w,ps(igrid)%x)
        end do
        mf_advance=.true.
        if(mype==0) then
          write(*,*) "itmf=",i
          write(*,*) '<CW sin theta>:',cwsin_theta_new
          write(*,*) '<f_i>:',f_i
          write(*,*) '----------------------------------------------------------'
        end if
      end if
      ! reconstruct AMR grid every 10 step
      if(mod(i,ditregrid)==0 .and. refine_max_level>1) call resettree
      if (i>=mf_it_max) then
        if(mod(i,10)/=0) then
          ! calculate metrics
          call metrics
          call printlog_mf
        end if
        if(mype==0) then
          write (*,*) 'Reach maximum iteration step!'
          write (*,*) 'The total iteration step is:', i
        end if
        exit
      end if
    enddo
    ! point bc mpi data type back to full type for MHD
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    bcphys=.true.
    ! set velocity back to zero and convert primitive variables back to conservative ones
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ps(igrid)%w(ixG^T,mom(:))=zero
       call phys_to_conserved(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
    end do
    global_time=tmpt
    it=tmpit
    if (mype==0) call MPI_FILE_CLOSE(fhmf,ierrmpi)
    mf_advance=.false.
    ! restore stagger_grid value
    stagger_grid=stagger_flag
    if(mype==0) write(*,*) 'Magnetofriction phase took : ',MPI_WTIME()-time_in,' sec'
    contains

      subroutine metrics

        sum_jbb_ipe = 0.d0
        sum_j_ipe = 0.d0
        sum_l_ipe = 0.d0
        f_i_ipe = 0.d0
        volumepe=0.d0
        dsurface=zero
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          block=>ps(igrid)
          dvolume(ixM^T)=block%dvolume(ixM^T)
          do idims=1,ndim
            hxM^LL=ixM^LL-kr(idims,^D);
            dsurface(ixM^T)=dsurface(ixM^T)+block%surfaceC(hxM^T,idims)
          end do
          ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
          call mask_inner(ixG^LL,ixM^LL,ps(igrid)%w,ps(igrid)%x)
          sum_jbb_ipe = sum_jbb_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
            ps(igrid)%x,1,patchwi)
          sum_j_ipe   = sum_j_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
            ps(igrid)%x,2,patchwi)
          f_i_ipe=f_i_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
            ps(igrid)%x,3,patchwi)
          sum_l_ipe   = sum_l_ipe+integral_grid_mf(ixG^LL,ixM^LL,ps(igrid)%w,&
            ps(igrid)%x,4,patchwi)
        end do
        call MPI_ALLREDUCE(sum_jbb_ipe,sum_jbb,1,MPI_DOUBLE_PRECISION,&
           MPI_SUM,icomm,ierrmpi)
        call MPI_ALLREDUCE(sum_j_ipe,sum_j,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        call MPI_ALLREDUCE(f_i_ipe,f_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        call MPI_ALLREDUCE(sum_l_ipe,sum_l,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        call MPI_ALLREDUCE(volumepe,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        ! current- and volume-weighted average of the sine of the angle
        ! between the magnetic field B and the current density J
        cwsin_theta_new = sum_jbb/sum_j
        ! volume-weighted average of the absolute value of the fractional
        ! magnetic flux change
        f_i = f_i/volume
        sum_j=sum_j/volume
        sum_l=sum_l/volume
      end subroutine metrics

      subroutine mask_inner(ixI^L,ixO^L,w,x)

        integer, intent(in)         :: ixI^L,ixO^L
        double precision, intent(in):: w(ixI^S,nw),x(ixI^S,1:ndim)
        double precision            :: xO^L
        integer                     :: ix^D

        {xOmin^D = xprobmin^D + 0.05d0*(xprobmax^D-xprobmin^D)\}
        {xOmax^D = xprobmax^D - 0.05d0*(xprobmax^D-xprobmin^D)\}
        if(slab) then
          xOmin^ND = xprobmin^ND
        else
          xOmin1 = xprobmin1
        end if

        {do ix^DB=ixOmin^DB,ixOmax^DB\}
            if({ x(ix^DD,^D) > xOmin^D .and. x(ix^DD,^D) < xOmax^D | .and. }) then
              patchwi(ix^D)=.true.
              volumepe=volumepe+dvolume(ix^D)
            else
              patchwi(ix^D)=.false.
            endif
        {end do\}

      end subroutine mask_inner

      subroutine printlog_mf
        integer :: amode, status(MPI_STATUS_SIZE)
        character(len=800) :: filename,filehead
        character(len=2048) :: line,datastr
        logical, save :: logmfopened=.false.

        if(mype==0) then
          if(.not.logmfopened) then
            ! generate filename
            write(filename,"(a,a)") TRIM(base_filename), "_mflog.csv"

            amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
            amode=ior(amode,MPI_MODE_APPEND)
            call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,fhmf,ierrmpi)
            logmfopened=.true.
            filehead="  itmf,  dt,  <f_i>,  <CW sin theta>,  <Current>,  <Lorenz force>"
            call MPI_FILE_WRITE(fhmf,filehead,len_trim(filehead), &
                                MPI_CHARACTER,status,ierrmpi)
            call MPI_FILE_WRITE(fhmf,achar(10),1,MPI_CHARACTER,status,ierrmpi)
          end if
          line=''
          write(datastr,'(i6,a)') i,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') dtfff,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') f_i,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') cwsin_theta_new,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') sum_j,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6)') sum_l
          line=trim(line)//trim(datastr)//new_line('A')
          call MPI_FILE_WRITE(fhmf,line,len_trim(line),MPI_CHARACTER,status,ierrmpi)
        end if

      end subroutine printlog_mf

      function integral_grid_mf(ixI^L,ixO^L,w,x,iw,patchwi)
        use mod_geometry

        integer, intent(in)                :: ixI^L,ixO^L,iw
        double precision, intent(in)       :: x(ixI^S,1:ndim)
        double precision, intent(in)       :: w(ixI^S,nw+nwauxio)
        logical, intent(in) :: patchwi(ixI^S)

        double precision, dimension(ixI^S,1:ndir) :: bvec,qvec,current
        double precision :: integral_grid_mf,tmp(ixI^S),b_mag(ixI^S)
        integer :: ix^D,i,idirmin,idir,jdir,kdir

        integral_grid_mf=0.d0
        select case(iw)
         case(1)
          ! Sum(dvolume*|JxB|/|B|)
          if(B0field) then
            bvec(ixI^S,:)=w(ixI^S,mag(:))+block%b0(ixI^S,mag(:),0)
          else
            bvec(ixI^S,:)=w(ixI^S,mag(:))
          endif
          call get_current(w,ixI^L,ixO^L,idirmin,current)
          ! calculate Lorentz force
          qvec(ixO^S,1:ndir)=zero
          do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
             if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
                if(lvc(idir,jdir,kdir)==1)then
                   qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
                else
                   qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
                endif
             endif
          enddo; enddo; enddo

          {do ix^DB=ixOmin^DB,ixOmax^DB\}
             if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+sqrt(sum(qvec(ix^D,:)**2)/&
                               sum(bvec(ix^D,:)**2))*dvolume(ix^D)
          {end do\}
         case(2)
          ! Sum(dvolume*|J|)
          call get_current(w,ixI^L,ixO^L,idirmin,current)
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
             if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+sqrt(sum(current(ix^D,:)**2))*&
                               dvolume(ix^D)
          {end do\}
         case(3)
          ! f_i solenoidal property of B: (dvolume |div B|)/(dsurface |B|)
          ! Sum(dvolume*f_i)
          if(B0field) then
            bvec(ixI^S,:)=w(ixI^S,mag(:))+block%b0(ixI^S,mag(:),0)
          else
            bvec(ixI^S,:)=w(ixI^S,mag(:))
          endif
          call divvector(bvec,ixI^L,ixO^L,tmp)
          {do ix^DB=ixOmin^DB,ixOmax^DB\}
             if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+abs(tmp(ix^D))*&
                   dvolume(ix^D)**2/sqrt(sum(bvec(ix^D,:)**2))/dsurface(ix^D)
          {end do\}
         case(4)
          ! Sum(|JxB|)
          if(B0field) then
            bvec(ixI^S,:)=w(ixI^S,mag(:))+block%b0(ixI^S,mag(:),0)
          else
            bvec(ixI^S,:)=w(ixI^S,mag(:))
          endif
          call curlvector(bvec,ixI^L,ixO^L,current,idirmin,1,ndir)
          ! calculate Lorentz force
          qvec(ixO^S,1:ndir)=zero
          do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
             if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixO^S)=current(ixO^S,jdir)*bvec(ixO^S,kdir)
                if(lvc(idir,jdir,kdir)==1)then
                   qvec(ixO^S,idir)=qvec(ixO^S,idir)+tmp(ixO^S)
                else
                   qvec(ixO^S,idir)=qvec(ixO^S,idir)-tmp(ixO^S)
                endif
             endif
          enddo; enddo; enddo

          {do ix^DB=ixOmin^DB,ixOmax^DB\}
             if(patchwi(ix^D)) integral_grid_mf=integral_grid_mf+sqrt(sum(qvec(ix^D,:)**2))*dvolume(ix^D)
          {end do\}
        end select
        return
      end function integral_grid_mf

  end subroutine magnetofriction

  subroutine mf_velocity_update(dtfff)
    use mod_global_parameters

    double precision, intent(in) :: dtfff
    integer :: i,iigrid, igrid
    double precision :: vhatmax,vhatmax_pe,vhatmaxgrid

    vhatmax_pe=smalldouble
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      call vhat(ps(igrid)%w,ps(igrid)%x,ixG^LL,ixM^LL,vhatmaxgrid)
      vhatmax_pe=max(vhatmax_pe,vhatmaxgrid)
    end do
    call MPI_ALLREDUCE(vhatmax_pe,vhatmax,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
                           icomm,ierrmpi)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      ! calculate frictional velocity
      call frictional_velocity(ps(igrid)%w,ps(igrid)%x,ixG^LL,ixM^LL,vhatmax,dtfff)
    end do

  end subroutine mf_velocity_update

  subroutine vhat(w,x,ixI^L,ixO^L,vhatmaxgrid)
    ! Calculate v_hat
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout)  :: w(ixI^S,nw)
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: vhatmaxgrid

    double precision :: current(ixI^S,7-2*ndir:3),tmp(ixI^S),dxhm
    double precision :: dxhms(ixO^S)
    integer :: idirmin,idir,jdir,kdir

    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixI^S,mom(:))=0.d0
    ! calculate Lorentz force
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
       if(lvc(idir,jdir,kdir)/=0)then
          if(B0field) then
            tmp(ixO^S)=current(ixO^S,jdir)*(w(ixO^S,mag(kdir))+block%b0(ixO^S,kdir,0))
          else
            tmp(ixO^S)=current(ixO^S,jdir)*w(ixO^S,mag(kdir))
          endif
          if(lvc(idir,jdir,kdir)==1)then
             w(ixO^S,mom(idir))=w(ixO^S,mom(idir))+tmp(ixO^S)
          else
             w(ixO^S,mom(idir))=w(ixO^S,mom(idir))-tmp(ixO^S)
          endif
       endif
    enddo; enddo; enddo

    if(B0field) then
      tmp(ixO^S)=sum((w(ixO^S,mag(:))+block%b0(ixO^S,:,0))**2,dim=ndim+1)  ! |B|**2
    else
      tmp(ixO^S)=sum(w(ixO^S,mag(:))**2,dim=ndim+1)         ! |B|**2
    endif

    if(slab_uniform) then
      dxhm=dble(ndim)/(^D&1.0d0/dxlevel(^D)+)
      do idir=1,ndir
        w(ixO^S,mom(idir))=dxhm*w(ixO^S,mom(idir))/tmp(ixO^S)
      end do
    else
      dxhms(ixO^S)=dble(ndim)/sum(1.d0/block%dx(ixO^S,:),dim=ndim+1)
      do idir=1,ndir
        w(ixO^S,mom(idir))=dxhms(ixO^S)*w(ixO^S,mom(idir))/tmp(ixO^S)
      end do
    end if
    vhatmaxgrid=maxval(sqrt(sum(w(ixO^S,mom(:))**2,dim=ndim+1)))

  end subroutine vhat

  subroutine frictional_velocity(w,x,ixI^L,ixO^L,qvmax,qdt)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim),qdt,qvmax
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: dxhm,disbd(6),bfzone^D
    double precision :: dxhms(ixO^S)
    integer :: ix^D, idir
    logical :: buffer

    if(slab_uniform) then
      dxhm=dble(ndim)/(^D&1.0d0/dxlevel(^D)+)
      dxhm=mf_cc*mf_cy/qvmax*dxhm/qdt
      ! dxhm=mf_cc*mf_cy/qvmax
      w(ixO^S,mom(:))=w(ixO^S,mom(:))*dxhm
    else
      dxhms(ixO^S)=dble(ndim)/sum(1.d0/block%dx(ixO^S,:),dim=ndim+1)
      dxhms(ixO^S)=mf_cc*mf_cy/qvmax*dxhms(ixO^S)/qdt
      ! dxhm=mf_cc*mf_cy/qvmax
      do idir=1,ndir
        w(ixO^S,mom(idir))=w(ixO^S,mom(idir))*dxhms(ixO^S)
      end do
    end if

{^IFTHREED
    bfzone1=0.05d0*(xprobmax1-xprobmin1)
    bfzone2=0.05d0*(xprobmax2-xprobmin2)
    bfzone3=0.05d0*(xprobmax3-xprobmin3)
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
       disbd(1)=x(ix^D,1)-xprobmin1
       disbd(2)=xprobmax1-x(ix^D,1)
       disbd(3)=x(ix^D,2)-xprobmin2
       disbd(4)=xprobmax2-x(ix^D,2)
       disbd(5)=x(ix^D,3)-xprobmin1
       disbd(6)=xprobmax3-x(ix^D,3)

       if(slab) then
         if(disbd(1)<bfzone1) then
           w(ix^D,mom(:))=(1.d0-((bfzone1-disbd(1))/bfzone1)**2)*w(ix^D,mom(:))
         endif
       else
         if(disbd(5)<bfzone3) then
           w(ix^D,mom(:))=(1.d0-((bfzone3-disbd(5))/bfzone3)**2)*w(ix^D,mom(:))
         endif
       end if
       if(disbd(2)<bfzone1) then
         w(ix^D,mom(:))=(1.d0-((bfzone1-disbd(2))/bfzone1)**2)*w(ix^D,mom(:))
       endif
       if(disbd(3)<bfzone2) then
         w(ix^D,mom(:))=(1.d0-((bfzone2-disbd(3))/bfzone2)**2)*w(ix^D,mom(:))
       endif
       if(disbd(4)<bfzone2) then
         w(ix^D,mom(:))=(1.d0-((bfzone2-disbd(4))/bfzone2)**2)*w(ix^D,mom(:))
       endif
       if(disbd(6)<bfzone3) then
         w(ix^D,mom(:))=(1.d0-((bfzone3-disbd(6))/bfzone3)**2)*w(ix^D,mom(:))
       endif
    {end do\}
}
  end subroutine frictional_velocity

  subroutine advectmf(idim^LIM,qt,qdt)
    !  integrate all grids by one step of its delta(global_time)
    ! This subroutine is in VAC terminology equivalent to
    ! `advect' (with the difference that it will `advect' all grids)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: idim^LIM
    double precision, intent(in) :: qt, qdt

    integer :: iigrid, igrid

    call init_comm_fix_conserve(idim^LIM,ndir)
    fix_conserve_at_step = mf_advance .and. levmax>levmin

    ! copy w instead of wold because of potential use of dimsplit or sourcesplit
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps1(igrid)%w=ps(igrid)%w
    end do

    istep=0

    select case (t_stepper)
     case (onestep)
       call advect1mf(flux_method,qdt,one,idim^LIM,qt,ps1,qt,ps)
     case (twostep)
       ! predictor step
       fix_conserve_at_step = .false.
       call advect1mf(typepred1,qdt,half,idim^LIM,qt,ps,qt,ps1)
       ! corrector step
       fix_conserve_at_step = mf_advance .and. levmax>levmin
       call advect1mf(flux_method,qdt,one,idim^LIM,qt+half*qdt,ps1,qt,ps)
     case (threestep)
       ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
       call advect1mf(flux_method,qdt,one,idim^LIM,qt,ps,qt,ps1)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixG^T,1:nwflux)=0.75d0*ps(igrid)%w(ixG^T,1:nwflux)+0.25d0*&
               ps1(igrid)%w(ixG^T,1:nwflux)
          if (nw>nwflux) ps2(igrid)%w(ixG^T,nwflux+1:nw) = &
               ps(igrid)%w(ixG^T,nwflux+1:nw)
       end do

       call advect1mf(flux_method,qdt,0.25d0,idim^LIM,qt+qdt,ps1,qt+dt*0.25d0,ps2)

       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ps(igrid)%w(ixG^T,1:nwflux)=1.0d0/3.0d0*ps(igrid)%w(ixG^T,1:nwflux)+&
               2.0d0/3.0d0*ps2(igrid)%w(ixG^T,1:nwflux)
       end do

       call advect1mf(flux_method,qdt,2.0d0/3.0d0,idim^LIM,qt+qdt/2.0d0,ps2,qt+qdt/3.0d0,ps)
     case default
       call mpistop("unkown time_stepper in advectmf")
    end select

  end subroutine advectmf

  subroutine advect1mf(method,dtin,dtfactor,idim^LIM,qtC,psa,qt,psb)
    ! Integrate all grids by one partial step
    ! This subroutine is equivalent to VAC's `advect1', but does
    ! the advection for all grids
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve

    integer, intent(in) :: idim^LIM
    type(state) :: psa(max_blocks)! Compute fluxes based on this state
    type(state) :: psb(max_blocks) ! update on this state
    double precision, intent(in) :: dtin,dtfactor, qtC, qt
    integer, intent(in) :: method(nlevelshi)

    double precision :: qdt
    integer :: iigrid, igrid, level, i^D
    logical :: setigrid

    istep=istep+1

    ! loop over all grids to arrive at equivalent
    qdt=dtfactor*dtin
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       block=>ps(igrid)
       level=node(plevel_,igrid)

       call process1_gridmf(method(level),igrid,qdt,ixG^LL,idim^LIM,qtC,&
                       psa(igrid)%w,qt,psb(igrid)%w,pso(igrid)%w)
    end do

    ! opedit: Send flux for all grids, expects sends for all
    ! nsend_fc(^D), set in connectivity.t.

    if (fix_conserve_at_step) then
      call recvflux(idim^LIM)
      call sendflux(idim^LIM)
      call fix_conserve(psb,idim^LIM,mag(1),ndir)
    end if
    ! point bc mpi datatype to partial type for magnetic field
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    ! update B in ghost cells
    call getbc(qt+qdt,qdt,psb,mag(1),ndir)
    ! calculate magnetofrictional velocity
    call mf_velocity_update(qdt)
    ! point bc mpi datatype to partial type for velocity field
    type_send_srl=>type_send_srl_p2
    type_recv_srl=>type_recv_srl_p2
    type_send_r=>type_send_r_p2
    type_recv_r=>type_recv_r_p2
    type_send_p=>type_send_p_p2
    type_recv_p=>type_recv_p_p2
    ! update magnetofrictional velocity in ghost cells
    call getbc(qt+qdt,qdt,psb,mom(1),ndir)

  end subroutine advect1mf

  subroutine process1_gridmf(method,igrid,qdt,ixG^L,idim^LIM,qtC,wCT,qt,w,wold)
    ! This subroutine is equivalent to VAC's `advect1' for one grid
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: method
    integer, intent(in) :: igrid, ixG^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision :: wCT(ixG^S,1:nw), w(ixG^S,1:nw), wold(ixG^S,1:nw)
    double precision :: dx^D, fC(ixG^S,1:ndir,1:ndim)
    integer :: ixO^L

    dx^D=rnode(rpdx^D_,igrid);
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    fC=0.d0

    ixO^L=ixG^L^LSUBnghostcells;
    select case (method)
     case (fs_cd4)
       !================================
       ! 4th order central difference
       !================================
       call centdiff4mf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,ps(igrid)%x)
     case (fs_tvdlf)
       !================================
       ! TVDLF
       !================================
       call tvdlfmf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,ps(igrid)%x)
     case (fs_hancock)
       ! hancock predict (first) step for twostep tvdlf and tvdmu scheme
       call hancockmf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,dx^D,ps(igrid)%x)
     case (fs_fd)
       !================================
       ! finite difference
       !================================
       call fdmf(qdt,ixG^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,ps(igrid)%x)
    case default
       call mpistop("unknown flux scheme in advect1_gridmf")
    end select

    if (fix_conserve_at_step) then
      call store_flux(igrid,fC,idim^LIM,ndir)
    end if

  end subroutine process1_gridmf

  subroutine upwindLRmf(ixI^L,ixL^L,ixR^L,idim,w,wCT,wLC,wRC,x)
    ! Determine the upwinded wLC(ixL) and wRC(ixR) from w.
    ! the wCT is only used when PPM is exploited.
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idim
    double precision, dimension(ixI^S,1:nw) :: w, wCT
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    double precision, dimension(ixI^S,1:ndim) :: x

    integer :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)

    if (type_limiter(block%level) == limiter_mp5) then
       call MP5limiter(ixI^L,ixL^L,idim,w,wLC,wRC)
    else if (type_limiter(block%level) == limiter_ppm) then
       call PPMlimiter(ixI^L,ixM^LL,idim,w,wCT,wLC,wRC)
    else
       jxR^L=ixR^L+kr(idim,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idim,^D);
       jxC^L=ixC^L+kr(idim,^D);

       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLC(ixL^S,iw)=dlog10(wLC(ixL^S,iw))
             wRC(ixR^S,iw)=dlog10(wRC(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idim,type_limiter(block%level),ldw,rdw)
          wLC(ixL^S,iw)=wLC(ixL^S,iw)+half*ldw(ixL^S)
          wRC(ixR^S,iw)=wRC(ixR^S,iw)-half*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
             wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
          end if
       end do

    endif

  end subroutine upwindLRmf

  subroutine getfluxmf(w,x,ixI^L,ixO^L,idir,idim,f)
    ! Calculate lux f_idim[idir] within ixO^L.
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idir, idim
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision,intent(out)    :: f(ixI^S)

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    if (idim==idir) then
      f(ixO^S)=zero
    else
      f(ixO^S)=w(ixO^S,mom(idim))*w(ixO^S,mag(idir))-w(ixO^S,mag(idim))*w(ixO^S,mom(idir))
      if (B0field) then
        f(ixO^S)=f(ixO^S)&
              +w(ixO^S,mom(idim))*block%B0(ixO^S,idir,idim)&
              -w(ixO^S,mom(idir))*block%B0(ixO^S,idim,idim)
      end if
    end if

  end subroutine getfluxmf

  subroutine tvdlfmf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,wold,fC,dx^D,x)
    ! method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
    ! method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.
    use mod_global_parameters

    double precision, intent(in)                         :: qdt, qtC, qt, dx^D
    integer, intent(in)                                  :: ixI^L, ixO^L, idim^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
    double precision, dimension(ixI^S,1:nw)               :: wCT, wnew, wold
    double precision, dimension(ixI^S,1:ndir,1:ndim)        :: fC

    double precision, dimension(ixI^S,1:nw) :: wLC, wRC, wmean
    double precision, dimension(ixI^S)      :: fLC, fRC
    double precision, dimension(ixI^S)      :: cmaxC
    double precision :: dxinv(1:ndim), inv_volume(ixO^S)
    integer :: idims, idir, ix^L, hxO^L, ixC^L, ixCR^L, jxC^L, kxC^L, kxR^L

    ! The flux calculation contracts by one in the idim direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
       call mpistop("Error in tvdlfmf: Nonconforming input limits")

    ^D&dxinv(^D)=-qdt/dx^D;
    fC=0.d0
    do idims= idim^LIM
       b0i=idims

       hxO^L=ixO^L-kr(idims,^D);
       ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2
       jxC^L=ixC^L+kr(idims,^D);
       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);
       ixCR^L=ixC^L;

       wRC(kxC^S,1:nwflux)=wCT(kxR^S,1:nwflux)
       wLC(kxC^S,1:nwflux)=wCT(kxC^S,1:nwflux)

       call upwindLRmf(ixI^L,ixCR^L,ixCR^L,idims,wCT,wCT,wLC,wRC,x)

       ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
       ! the maximum eigenvalue, it is calculated in advance.
       ! determine mean state and store in wLC
       wmean=0.5d0*(wLC+wRC)
       call getcmaxfff(wmean,ixG^LL,ixC^L,idims,cmaxC)

       ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each idir
       do idir=1,ndir
          call getfluxmf(wLC,x,ixG^LL,ixC^L,idir,idims,fLC)
          call getfluxmf(wRC,x,ixG^LL,ixC^L,idir,idims,fRC)
          ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
          fLC(ixC^S)=half*(fLC(ixC^S)+fRC(ixC^S))

          ! Add TVDLF dissipation to the flux
          if (idir==idims) then
            fLC(ixC^S)=fLC(ixC^S)-mf_tvdlfeps*tvdlfeps*cmaxC(ixC^S)*half*(wRC(ixC^S,mag(idir))-wLC(ixC^S,mag(idir)))
          end if
          if (slab_uniform) then
            fC(ixC^S,idir,idims)=fLC(ixC^S)
          else
            fC(ixC^S,idir,idims)=block%surfaceC(ixC^S,idims)*fLC(ixC^S)
          end if

       end do ! Next idir
    end do ! Next idims
    b0i=0

    !Now update the state:
    do idims= idim^LIM
       hxO^L=ixO^L-kr(idims,^D);
       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          fC(ixI^S,:,idims)=dxinv(idims)*fC(ixI^S,:,idims)
          wnew(ixO^S,mag(:))=wnew(ixO^S,mag(:)) &
               + (fC(ixO^S,:,idims)-fC(hxO^S,:,idims))
       else
          inv_volume = 1.0d0/block%dvolume(ixO^S)
          fC(ixI^S,:,idims)=-qdt*fC(ixI^S,:,idims)

          do idir = 1, ndir
             wnew(ixO^S,mag(idir))=wnew(ixO^S,mag(idir)) + (fC(ixO^S,idir,idims)-fC(hxO^S,idir,idims)) * &
                  inv_volume
          end do
       end if

    end do ! Next idims

    if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call divbclean(qdt,ixI^L,ixO^L,wCT,wnew,x)

  end subroutine tvdlfmf

  subroutine hancockmf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,dx^D,x)
    ! The non-conservative Hancock predictor for TVDLFmf
    ! on entry:
    ! input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
    ! one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
    ! on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), wnew(ixI^S,1:nw)

    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    double precision, dimension(ixI^S) :: fLC, fRC
    double precision :: dxinv(1:ndim)
    integer :: idims, idir, ix^L, hxO^L, ixtest^L

    ! Expand limits in each idims direction in which fluxes are added
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADDkr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
       call mpistop("Error in Hancockmf: Nonconforming input limits")

    ^D&dxinv(^D)=-qdt/dx^D;
    do idims= idim^LIM
       b0i=idims
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxO^L=ixO^L-kr(idims,^D);

       wRC(hxO^S,1:nwflux)=wCT(ixO^S,1:nwflux)
       wLC(ixO^S,1:nwflux)=wCT(ixO^S,1:nwflux)

       call upwindLRmf(ixI^L,ixO^L,hxO^L,idims,wCT,wCT,wLC,wRC,x)

       ! Advect mag(idir)
       do idir=1,ndir
          ! Calculate the fLC and fRC fluxes
          call getfluxmf(wRC,x,ixI^L,hxO^L,idir,idims,fRC)
          call getfluxmf(wLC,x,ixI^L,ixO^L,idir,idims,fLC)

          if (slab_uniform) then
             wnew(ixO^S,mag(idir))=wnew(ixO^S,mag(idir))+dxinv(idims)* &
                              (fLC(ixO^S)-fRC(hxO^S))
          else
             wnew(ixO^S,mag(idir))=wnew(ixO^S,mag(idir))-qdt/block%dvolume(ixO^S) &
                   *(block%surfaceC(ixO^S,idims)*fLC(ixO^S) &
                    -block%surfaceC(hxO^S,idims)*fRC(hxO^S))
          end if
       end do
    end do ! next idims
    b0i=0

    if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,wnew,x)

  end subroutine hancockmf

  subroutine fdmf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,wnew,wold,fC,dx^D,x)
    use mod_global_parameters
    double precision, intent(in)                                     :: qdt, qtC, qt, dx^D
    integer, intent(in)                                              :: ixI^L, ixO^L, idim^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in)            :: x

    double precision, dimension(ixI^S,1:nw), intent(inout)           :: wCT, wnew, wold
    double precision, dimension(ixI^S,1:ndir,1:ndim), intent(out)  :: fC

    double precision, dimension(ixI^S)                               :: fCT
    double precision, dimension(ixI^S,1:nw)                          :: fm, fp, fmR, fpL
    double precision, dimension(ixI^S)                               :: v
    double precision                                                 :: dxinv(1:ndim)
    integer                                                          :: idims, idir, ixC^L, ix^L, hxO^L, ixCR^L

    ^D&dxinv(^D)=-qdt/dx^D;
    do idims= idim^LIM

       ! Get fluxes for the whole grid (mesh+nghostcells)
       {^D& ixCmin^D = ixOmin^D - nghostcells * kr(idims,^D)\}
       {^D& ixCmax^D = ixOmax^D + nghostcells * kr(idims,^D)\}

       hxO^L=ixO^L-kr(idims,^D);
       ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixmax^D=ixOmax^D; ixmin^D=hxOmin^D;
       ixCR^L=ixC^L;

       do idir=1,ndir
          call getfluxmf(wCT,x,ixG^LL,ixCR^L,idir,idims,fCT)
          ! Lax-Friedrich splitting:
          fp(ixCR^S,mag(idir)) = half * (fCT(ixCR^S) + mf_tvdlfeps * tvdlfeps * cmax_global * wCT(ixCR^S,mag(idir)))
          fm(ixCR^S,mag(idir)) = half * (fCT(ixCR^S) - mf_tvdlfeps * tvdlfeps * cmax_global * wCT(ixCR^S,mag(idir)))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructLmf(ixI^L,ix^L,idims,fp,fpL)
       call reconstructRmf(ixI^L,ix^L,idims,fm,fmR)

       do idir=1,ndir
          if (slab_uniform) then
             fC(ix^S,idir,idims) = dxinv(idims) * (fpL(ix^S,mag(idir)) + fmR(ix^S,mag(idir)))
             wnew(ixO^S,mag(idir))=wnew(ixO^S,mag(idir))+ &
                  (fC(ixO^S,idir,idims)-fC(hxO^S,idir,idims))
          else
             fC(ix^S,idir,idims)=-qdt*block%surfaceC(ix^S,idims) * (fpL(ix^S,mag(idir)) + fmR(ix^S,mag(idir)))
             wnew(ixO^S,mag(idir))=wnew(ixO^S,mag(idir))+ &
               (fC(ixO^S,idir,idims)-fC(hxO^S,idir,idims))/block%dvolume(ixO^S)
          end if
       end do ! iw loop

    end do !idims loop

    if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call divbclean(qdt,ixI^L,ixO^L,wCT,wnew,x)

  end subroutine fdmf

  subroutine reconstructLmf(ixI^L,iL^L,idims,w,wLC)
    use mod_global_parameters
    use mod_limiter

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(out)   :: wLC(ixI^S,1:nw)

    double precision                :: ldw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, iw

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterL(ixI^L,iL^L,idims,w,wLC)
    case (limiter_weno5)
       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,1)
    case (limiter_wenoz5)
       call WENO5limiterL(ixI^L,iL^L,idims,w,wLC,2)
    case default

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);

       wLC(kxC^S,1:nwflux) = w(kxC^S,1:nwflux)

       jxR^L=iL^L+kr(idims,^D);

       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,nwflux
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),ldw=ldw)

          wLC(iL^S,iw)=wLC(iL^S,iw)+half*ldw(iL^S)
       end do
    end select

  end subroutine reconstructLmf

  subroutine reconstructRmf(ixI^L,iL^L,idims,w,wRC)
    use mod_global_parameters
    use mod_limiter

    integer, intent(in)             :: ixI^L, iL^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision, intent(out)   :: wRC(ixI^S,1:nw)

    double precision                :: rdw(ixI^S), dwC(ixI^S)
    integer                         :: jxR^L, ixC^L, jxC^L, kxC^L, kxR^L, iw

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterR(ixI^L,iL^L,idims,w,wRC)
    case (limiter_weno5)
       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,1)
    case (limiter_wenoz5)
       call WENO5limiterR(ixI^L,iL^L,idims,w,wRC,2)
    case default

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       wRC(kxC^S,1:nwflux)=w(kxR^S,1:nwflux)

       jxR^L=iL^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=iLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,nwflux
          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,type_limiter(block%level),rdw=rdw)

          wRC(iL^S,iw)=wRC(iL^S,iw)-half*rdw(jxR^S)
       end do
    end select

  end subroutine reconstructRmf

  subroutine centdiff4mf(qdt,ixI^L,ixO^L,idim^LIM,qtC,wCT,qt,w,wold,fC,dx^D,x)
    ! Advance the flow variables from global_time to global_time+qdt within ixO^L by
    ! fourth order centered differencing in space
    ! for the dw/dt+dF_i(w)/dx_i=S type equation.
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idim^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D
    double precision :: wCT(ixI^S,1:nw), w(ixI^S,1:nw), wold(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: fC(ixI^S,1:ndir,1:ndim)

    double precision :: v(ixI^S,ndim), f(ixI^S)
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    double precision, dimension(ixI^S)      :: vLC, vRC,cmaxLC,cmaxRC
    double precision :: dxinv(1:ndim)
    integer :: idims, idir, idirmin,ix^D
    integer :: ix^L, hxO^L, ixC^L, jxC^L, hxC^L, kxC^L, kkxC^L, kkxR^L

    ! two extra layers are needed in each direction for which fluxes are added.
    ix^L=ixO^L;
    do idims= idim^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do

    if (ixI^L^LTix^L|.or.|.or.) then
       call mpistop("Error in evolve_CentDiff4: Non-conforming input limits")
    end if
    ^D&dxinv(^D)=-qdt/dx^D;

    ! Add fluxes to w
    do idims= idim^LIM
       b0i=idims
       ix^L=ixO^L^LADD2*kr(idims,^D);
       hxO^L=ixO^L-kr(idims,^D);

       ixCmin^D=hxOmin^D; ixCmax^D=ixOmax^D;
       hxC^L=ixC^L-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);
       kxC^L=ixC^L+2*kr(idims,^D);

       kkxCmin^D=ixImin^D; kkxCmax^D=ixImax^D-kr(idims,^D);
       kkxR^L=kkxC^L+kr(idims,^D);
       wRC(kkxC^S,1:nwflux)=wCT(kkxR^S,1:nwflux)
       wLC(kkxC^S,1:nwflux)=wCT(kkxC^S,1:nwflux)

       call upwindLRmf(ixI^L,ixC^L,ixC^L,idims,wCT,wCT,wLC,wRC,x)

       ! Calculate velocities from upwinded values
       call getcmaxfff(wLC,ixG^LL,ixC^L,idims,cmaxLC)
       call getcmaxfff(wRC,ixG^LL,ixC^L,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixC^S)=max(cmaxRC(ixC^S),cmaxLC(ixC^S))

       do idir=1,ndir
          ! Get non-transported flux
          call getfluxmf(wCT,x,ixI^L,ix^L,idir,idims,f)
          ! Center flux to interface
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixC^S,idir,idims)=(-f(kxC^S)+7.0d0*(f(jxC^S)+f(ixC^S))-f(hxC^S))/12.0d0
          ! add rempel dissipative flux, only second order version for now
          ! one could gradually reduce the dissipative flux to improve solutions
          ! for computing steady states (Keppens et al. 2012, JCP)
          fC(ixC^S,idir,idims)=fC(ixC^S,idir,idims)-mf_tvdlfeps*tvdlfeps*half*vLC(ixC^S) &
                                         *(wRC(ixC^S,mag(idir))-wLC(ixC^S,mag(idir)))

          if (slab_uniform) then
             fC(ixC^S,idir,idims)=dxinv(idims)*fC(ixC^S,idir,idims)
             ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)+f_(i-1))-f_(i-2)]/12
             w(ixO^S,mag(idir))=w(ixO^S,mag(idir))+(fC(ixO^S,idir,idims)-fC(hxO^S,idir,idims))
          else
             fC(ixC^S,idir,idims)=-qdt*block%surfaceC(ixC^S,idims)*fC(ixC^S,idir,idims)
             w(ixO^S,mag(idir))=w(ixO^S,mag(idir))+ &
                  (fC(ixO^S,idir,idims)-fC(hxO^S,idir,idims))/block%dvolume(ixO^S)
          end if
       end do    !next idir
    end do       !next idims
    b0i=0

    if (.not.slab) call addgeometrymf(qdt,ixI^L,ixO^L,wCT,w,x)
    call divbclean(qdt,ixI^L,ixO^L,wCT,w,x)

  end subroutine centdiff4mf

  subroutine getdtfff_courant(w,x,ixI^L,ixO^L,dtnew)
    ! compute CFL limited dt (for variable time stepping)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw), dtnew

    double precision :: courantmax, dxinv(1:ndim)
    double precision :: cmax(ixI^S),tmp(ixI^S),alfven(ixI^S)
    integer :: idims

    dtnew=bigdouble
    courantmax=0.d0
    ^D&dxinv(^D)=1.d0/dxlevel(^D);

    do idims=1,ndim
       call getcmaxfff(w,ixI^L,ixO^L,idims,cmax)
       cmax_mype = max(cmax_mype,maxval(cmax(ixO^S)))
       if (.not.slab_uniform) then
          tmp(ixO^S)=cmax(ixO^S)/block%dx(ixO^S,idims)
          courantmax=max(courantmax,maxval(tmp(ixO^S)))
       else
          tmp(ixO^S)=cmax(ixO^S)*dxinv(idims)
          courantmax=max(courantmax,maxval(tmp(ixO^S)))
       end if
    end do
    ! courantmax='max( c/dx)'
    if (courantmax>smalldouble)  dtnew=min(dtnew,mf_cc/courantmax)

  end subroutine getdtfff_courant

  subroutine getcmaxfff(w,ixI^L,ixO^L,idims,cmax)
    use mod_global_parameters

    logical :: new_cmax,needcmin
    integer, intent(in) :: ixI^L, ixO^L, idims
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(out) :: cmax(ixI^S)

    ! calculate alfven speed
    if(B0field) then
      cmax(ixO^S)=sqrt(sum((w(ixO^S,mag(:))+block%b0(ixO^S,:,0))**2,dim=ndim+1)/w(ixO^S,rho_))
    else
      cmax(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2,dim=ndim+1)/w(ixO^S,rho_))
    endif
    cmax(ixO^S)=cmax(ixO^S)+abs(w(ixO^S,mom(idims)))

  end subroutine getcmaxfff

  !> Clean divergence of magnetic field by Janhunen's and Linde's source terms
  subroutine divbclean(qdt,ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim),wCT(ixI^S,1:nw),qdt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: idims, ix^L, ixp^L, i^D, iside
    double precision :: divb(ixI^S),graddivb(ixI^S),bdivb(ixI^S,1:ndir)

    ! Calculate div B
    ix^L=ixO^L^LADD1;
    call get_divb(wCT,ixI^L,ix^L,divb)

    ixp^L=ixO^L;

    ! Add Linde's diffusive terms
    do idims=1,ndim
       ! Calculate grad_idim(divb)
       call gradient(divb,ixI^L,ixp^L,idims,graddivb)

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixp^S)=graddivb(ixp^S)*mf_cdivb/(^D&1.0d0/dxlevel(^D)**2+)
       else
          graddivb(ixp^S)=graddivb(ixp^S)*mf_cdivb &
                          /(^D&1.0d0/block%dx(ixp^S,^D)**2+)
       end if
       ! B_idim += eta*grad_idim(divb)
       ! with Janhunen's term
       !w(ixp^S,mag(idims))=w(ixp^S,mag(idims))+&
       !      graddivb(ixp^S)-qdt*w(ixp^S,mom(idims))*divb(ixp^S)
       ! without Janjunen's term
       w(ixp^S,mag(idims))=w(ixp^S,mag(idims))+&
             graddivb(ixp^S)
    end do

  end subroutine divbclean

  subroutine addgeometrymf(qdt,ixI^L,ixO^L,wCT,w,x)
    ! Add geometrical source terms to w
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)
    !.. local ..
    double precision :: tmp(ixI^S)
    integer          :: iw
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (coordinate)
    case (cylindrical)
      if(phi_>0) then
        ! s[Bphi]=(Bphi*vr-Br*vphi)/radius
        tmp(ixO^S)=(wCT(ixO^S,bphi_)*wCT(ixO^S,mom(1)) &
                   -wCT(ixO^S,br_)*wCT(ixO^S,mom(3)))
        w(ixO^S,bphi_)=w(ixO^S,bphi_)+qdt*tmp(ixO^S)/x(ixO^S,1)
      end if
    case (spherical)
    {^NOONED
      ! s[b2]=(vr*Btheta-vtheta*Br)/r
      !       + cot(theta)*psi/r
      tmp(ixO^S)= wCT(ixO^S,mom(1))*wCT(ixO^S,mag(2)) &
                 -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(1))
      if (B0field) then
         tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(1))*block%b0(ixO^S,2,0) &
                    -wCT(ixO^S,mom(2))*block%b0(ixO^S,1,0)
      end if
      ! Divide by radius and add to w
      w(ixO^S,mag(2))=w(ixO^S,mag(2))+qdt*tmp(ixO^S)/x(ixO^S,1)
    }
      if(ndir==3) then
        ! s[b3]=(vr*Bphi-vphi*Br)/r
        ! -cot(theta)*(vphi*Btheta-vtheta*Bphi)/r
        tmp(ixO^S)=wCT(ixO^S,mom(1))*wCT(ixO^S,mag(3)) &
                -wCT(ixO^S,mom(3))*wCT(ixO^S,mag(1)){^NOONED &
               -(wCT(ixO^S,mom(3))*wCT(ixO^S,mag(2)) &
                -wCT(ixO^S,mom(2))*wCT(ixO^S,mag(3)))*dcos(x(ixO^S,2)) &
                              /dsin(x(ixO^S,2)) }
        if (B0field) then
           tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,mom(1))*block%b0(ixO^S,3,0) &
              -wCT(ixO^S,mom(3))*block%b0(ixO^S,1,0){^NOONED &
              -(wCT(ixO^S,mom(3))*block%b0(ixO^S,2,0) &
               -wCT(ixO^S,mom(2))*block%b0(ixO^S,3,0))*dcos(x(ixO^S,2)) &
                              /dsin(x(ixO^S,2)) }
        end if
        ! Divide by radius and add to w
        w(ixO^S,mag(3))=w(ixO^S,mag(3))+qdt*tmp(ixO^S)/x(ixO^S,1)
      end if
    end select

  end subroutine addgeometrymf

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixI^L,ixO^L,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer :: idirmin0
    integer :: ixO^L, idirmin, ixI^L
    double precision :: w(ixI^S,1:nw)
    integer :: idir

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixI^S,7-2*ndir:3),bvec(ixI^S,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))

    call curlvector(bvec,ixI^L,ixO^L,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixO^S,idirmin0:3)=current(ixO^S,idirmin0:3)+&
        block%J0(ixO^S,idirmin0:3)

  end subroutine get_current

  !> Calculate div B within ixO
  subroutine get_divb(w,ixI^L,ixO^L,divb)
    use mod_geometry
    use mod_global_parameters

    integer, intent(in)                :: ixI^L, ixO^L
    double precision, intent(in)       :: w(ixI^S,1:nw)
    double precision                   :: divb(ixI^S)

    double precision                   :: bvec(ixI^S,1:ndir)

    bvec(ixI^S,:)=w(ixI^S,mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixI^L,ixO^L,divb)
    case("limited")
      call divvectorS(bvec,ixI^L,ixO^L,divb)
    end select
  end subroutine get_divb

end module mod_magnetofriction
