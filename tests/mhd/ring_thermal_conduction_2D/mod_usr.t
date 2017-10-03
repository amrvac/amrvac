! setup.pl -d=2
! test thermal conduction in a ring
module mod_usr
  use mod_mhd
  implicit none
  double precision :: integral_time=0.d0, k_perp

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_init_one_grid   => initonegrid_usr
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_special_convert => userspecialconvert
    !usr_process_global  => usrprocess_global

    call set_coordinate_system("Cartesian")
    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: r(ixI^S), theta(ixI^S), B(ixI^S)

    ! set all velocity to zero
    w(ixO^S, mom(:)) = zero
    ! uniform pressure
    w(ixO^S,rho_) =1.d0
    r(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where(x(ixO^S,1)>0.d0)
      theta(ixO^S)=atan(x(ixO^S,2)/x(ixO^S,1))
    elsewhere(x(ixO^S,1)<0.d0)
      theta(ixO^S)=dpi-atan(x(ixO^S,2)/abs(x(ixO^S,1)))
    elsewhere
      theta(ixO^S)=0.d0
    endwhere
    ! hot central circular spot with uniform pressure
    where(r(ixO^S)>0.5d0 .and. r(ixO^S)<0.7d0 .and. theta(ixO^S)>11.d0/12.d0*dpi .and. theta(ixO^S)<13.d0/12.d0*dpi)
      w(ixO^S,p_)=12.d0
    elsewhere
      w(ixO^S,p_)=10.d0
    endwhere
    ! straight line current
    B(ixO^S)=Busr/r(ixO^S)
    w(ixO^S,mag(1))=B(ixO^S)*dcos(theta(ixO^S)+0.5*dpi)
    w(ixO^S,mag(2))=B(ixO^S)*dsin(theta(ixO^S)+0.5*dpi)

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: r(ixI^S), pth(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    pth(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)

    r(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where(r(ixO^S)>0.5d0 .and. r(ixO^S)<0.7d0)
      w(ixO^S,nw+1)=dabs(pth(ixO^S)-61.d0/6.d0)
    elsewhere
      w(ixO^S,nw+1)=dabs(pth(ixO^S)-10.d0)
    endwhere
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='L'

  end subroutine specialvarnames_output

  subroutine userspecialconvert(qunitconvert)
    use mod_global_parameters
    integer, intent(in) :: qunitconvert

    call spatial_integral

  end subroutine userspecialconvert

  subroutine spatial_integral
    use mod_global_parameters

    double precision :: dvolume(ixG^T), timephy,xmom,ymom
    double precision, allocatable :: integral_ipe(:), integral_w(:)

    integer           :: nregions,ireg,i
    integer           :: iigrid,igrid,status(MPI_STATUS_SIZE),ni
    character(len=100):: filename,region
    character(len=1024) :: line, datastr
    logical           :: patchwi(ixG^T),alive

    ! number of integrals to perform
    ni=5
    allocate(integral_ipe(ni),integral_w(ni))
    integral_ipe=0.d0
    integral_w=0.d0
    integral_ipe(5)=bigdouble
    integral_w(5)=bigdouble
    
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>pw(igrid)
      if(slab) then
        dvolume(ixM^T)={rnode(rpdx^D_,igrid)|*}
      else
        dvolume(ixM^T)=block%dvolume(ixM^T)
      end if
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
      patchwi(ixM^T)=.true.
      integral_ipe(1)=integral_ipe(1)+ &
                integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,pw(igrid)%x,dvolume,1,patchwi)
      integral_ipe(2)=integral_ipe(2)+ &
                integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,pw(igrid)%x,dvolume,2,patchwi)
      integral_ipe(3)=max(integral_ipe(3), &
                integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,pw(igrid)%x,dvolume,3,patchwi))
      integral_ipe(4)=max(integral_ipe(4), &
                integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,pw(igrid)%x,dvolume,4,patchwi))
      integral_ipe(5)=min(integral_ipe(5), &
                integral_grid(ixG^LL,ixM^LL,pw(igrid)%w,pw(igrid)%x,dvolume,5,patchwi))
    end do
    call MPI_ALLREDUCE(integral_ipe(1:2),integral_w(1:2),2,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
    call MPI_ALLREDUCE(integral_ipe(3:4),integral_w(3:4),2,MPI_DOUBLE_PRECISION,&
                         MPI_MAX,icomm,ierrmpi)
    call MPI_ALLREDUCE(integral_ipe(5),integral_w(5),1,MPI_DOUBLE_PRECISION,&
                         MPI_MIN,icomm,ierrmpi)
    integral_w(1:2)=integral_w(1:2)/(^D&(xprobmax^D-xprobmin^D)*)
    integral_w(2)=dsqrt(integral_w(2))
    timephy=global_time
    if(mype==0) then
      write(filename,"(a,a)") TRIM(base_filename),"Lerrors.csv"
      inquire(file=filename,exist=alive)
      if(alive) then
        open(unit=21,file=filename,form='formatted',status='old',access='append')
      else
        open(unit=21,file=filename,form='formatted',status='new')
        write(21,'(a)') 'time, L1 error, L2 error, L infinity error, Tmax, Tmin' 
      endif
      write(datastr,'(es11.4, 2a)') timephy,', '
      line=datastr
      do i=1,ni-1
        write(datastr,"(es12.5, 2a)") integral_w(i),', '
        line = trim(line)//trim(datastr)
      end do
      write(datastr,"(es12.5)") integral_w(ni)
      line = trim(line)//trim(datastr)
      write(21,'(a)') trim(line)
      close(21)
    endif
    
    deallocate(integral_ipe,integral_w)
  
  end subroutine spatial_integral

  double precision function integral_grid(ixI^L,ixO^L,w,x,dvolume,intval,patchwi)
    use mod_global_parameters

    integer, intent(in)                :: ixI^L,ixO^L,intval
    double precision, intent(in)       :: x(ixI^S,1:ndim),dvolume(ixI^S)
    double precision, intent(in)       :: w(ixI^S,nw)
    logical, intent(in) :: patchwi(ixI^S)

    double precision :: r(ixI^S), pth(ixI^S),L(ixI^S)
    integer :: ix^D

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    pth(ixO^S)=pth(ixO^S)/w(ixO^S,rho_)

    r(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where(r(ixO^S)>0.5d0 .and. r(ixO^S)<0.7d0)
      L(ixO^S)=dabs(pth(ixO^S)-61.d0/6.d0)
    elsewhere
      L(ixO^S)=dabs(pth(ixO^S)-10.d0)
    endwhere


    integral_grid=0.d0
    select case(intval)
     case(1)
      ! L1
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+L(ix^D)*dvolume(ix^D)
      {end do\}
     case(2)
      ! L2
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
         if(patchwi(ix^D)) integral_grid=integral_grid+L(ix^D)**2*dvolume(ix^D)
      {end do\}
     case(3)
      ! L infnity
      integral_grid=maxval(L(ixO^S))
     case(4)
      ! T max
      integral_grid=maxval(pth(ixO^S))
     case(5)
      ! T min 
      integral_grid=minval(pth(ixO^S))
     case default
         call mpistop("intval not defined")
    end select

  end function integral_grid

  subroutine usrprocess_global(iit,qt)
    use mod_global_parameters
    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    double precision :: integral_g(2), integral_ipe(2)
    integer:: iigrid, igrid

    integral_ipe=0.d0
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
       block=>pw(igrid)
       call usrprocess_grid(ixG^LL,ixM^LL,qt,pw(igrid)%w,pw(igrid)%x,integral_g)
       integral_ipe=integral_ipe+integral_g
    end do
    call MPI_ALLREDUCE(integral_ipe,integral_g,2,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,icomm,ierrmpi)
    integral_time=integral_time+dt*integral_g(1)

    k_perp=integral_g(2)/integral_time

    if(timetosaveusr(2).and.mype==0) write(*,*) "K_perp", k_perp, k_perp/0.01d0, qt

  end subroutine usrprocess_global

  logical function timetosaveusr(ifile)
    use mod_global_parameters

    integer:: ifile
    logical:: oksave

    oksave=.false.
    if (global_time>=tsavelast(ifile)+dtsave(ifile)-smalldouble)then
       oksave=.true.
    endif

    timetosaveusr=oksave

    return
  end function timetosaveusr

  subroutine usrprocess_grid(ixI^L,ixO^L,qt,w,x,integral)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision, dimension(ixI^S,1:ndim) :: qvec, gradT
    double precision, dimension(ixI^S) :: tmp1, Te, minq, maxq, qd, r, theta
    double precision :: dxinv(ndim), alpha, integral(2)
    integer :: ix^D, ix^L, ixC^L, ixA^L, ixB^L, idims

    alpha=0.75d0
    ix^L=ixO^L^LADD1;

    dxinv=1.d0/dxlevel

    call mhd_get_pthermal(w,x,ixI^L,ixI^L,tmp1)
    ! compute the temperature
    Te(ixI^S)=tmp1(ixI^S)/w(ixI^S,rho_)
    ! ixC is cell-corner index
    ixCmax^D=ixOmax^D; ixCmin^D=ixOmin^D-1;
    ! T gradient at cell faces
    gradT=0.d0
    do idims=1,ndim
      ixBmin^D=ixmin^D;
      ixBmax^D=ixmax^D-kr(idims,^D);
      ixA^L=ixB^L+kr(idims,^D);
      gradT(ixB^S,idims)=(Te(ixA^S)-Te(ixB^S))*dxinv(idims)
    end do
    ! calculate thermal conduction flux with slope-limited symmetric scheme
    qvec=0.d0
    do idims=1,ndim
      ixB^L=ixO^L-kr(idims,^D);
      ixAmax^D=ixOmax^D; ixAmin^D=ixBmin^D;
      ! limited normal component
      minq(ixA^S)=min(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
      maxq(ixA^S)=max(alpha*gradT(ixA^S,idims),gradT(ixA^S,idims)/alpha)
      ! eq (19)
      qd=0.d0
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixCmin^D+ix^D;
           ixBmax^D=ixCmax^D+ix^D;
           qd(ixC^S)=qd(ixC^S)+gradT(ixB^S,idims)
         end if
      {end do\}
      ! temperature gradient at cell corner
      qd(ixC^S)=qd(ixC^S)*0.5d0**(ndim-1)
      ! eq (21)
      {do ix^DB=0,1 \}
         if({ ix^D==0 .and. ^D==idims | .or.}) then
           ixBmin^D=ixAmin^D-ix^D;
           ixBmax^D=ixAmax^D-ix^D;
           where(qd(ixB^S)<=minq(ixA^S))
             qd(ixB^S)=minq(ixA^S)
           elsewhere(qd(ixB^S)>=maxq(ixA^S))
             qd(ixB^S)=maxq(ixA^S)
           end where
           qvec(ixA^S,idims)=qvec(ixA^S,idims)+qd(ixB^S)
         end if
      {end do\}
      qvec(ixA^S,idims)=qvec(ixA^S,idims)*0.5d0**(ndim-1)
    end do

    qd=0.d0
    do idims=1,ndim
      qvec(ix^S,idims)=dxinv(idims)*qvec(ix^S,idims)
      ixB^L=ixO^L-kr(idims,^D);
      qd(ixO^S)=qd(ixO^S)+qvec(ixO^S,idims)-qvec(ixB^S,idims)
    end do

    r(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    where(x(ixO^S,1)>0.d0)
      theta(ixO^S)=atan(x(ixO^S,2)/x(ixO^S,1))
    elsewhere(x(ixO^S,1)<0.d0)
      theta(ixO^S)=dpi-atan(x(ixO^S,2)/abs(x(ixO^S,1)))
    elsewhere
      theta(ixO^S)=0.d0
    endwhere
    alpha=product(dxlevel(:))
    integral=0.d0
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
      if(r(ix^D)>0.5d0 .and. r(ix^D)<0.7d0) then
        integral(1)=integral(1)+qd(ix^D)*alpha
        if(theta(ix^D)>11.d0/12.d0*dpi .and. theta(ix^D)<13.d0/12.d0*dpi) then 
          integral(2)=integral(2)+(Te(ix^D)-12.d0)*alpha
        else
          integral(2)=integral(2)+(Te(ix^D)-10.d0)*alpha
        end if
      end if
    {end do\}

  end subroutine usrprocess_grid

end module mod_usr
