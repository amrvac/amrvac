module mod_errest
  use mod_comm_lib, only: mpistop
  implicit none
  private

  public :: errest, forcedrefine_grid_io 
 
contains

  !> Do all local error estimation which determines (de)refinement
  subroutine errest
    use mod_forest, only: refine, buffer
    use mod_global_parameters

    integer :: igrid, iigrid, ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,&
       ixCoGmax2,ixCoGmax3
    double precision :: factor
    logical, dimension(:,:), allocatable :: refine2

    if (igridstail==0) return

    select case (refine_criterion)
    case (1) 
       ! all refinement solely based on user routine usr_refine_grid
    case (2)
       ! Error estimation is based on Lohner's original scheme
    !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          block=>ps(igrid)
          call lohner_orig_grid(igrid)
       end do
    !$OMP END PARALLEL DO
    case (3)
       ! Error estimation is based on Lohner's scheme
    !$OMP PARALLEL DO PRIVATE(igrid)
       do iigrid=1,igridstail; igrid=igrids(iigrid);
          block=>ps(igrid)
          call lohner_grid(igrid)
       end do
    !$OMP END PARALLEL DO
    case default
       call mpistop("Unknown error estimator")
    end select

    ! enforce additional refinement on e.g. coordinate and/or time info here
    if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
       allocate(refine2(max_blocks,npe))
       call MPI_ALLREDUCE(refine,refine2,max_blocks*npe,MPI_LOGICAL,MPI_LOR,&
           icomm,ierrmpi)
       refine=refine2
    end if
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       block=>ps(igrid)
       call forcedrefine_grid(igrid,ps(igrid)%w)
    end do
    !$OMP END PARALLEL DO

    if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) buffer=.false.

  end subroutine errest

  subroutine lohner_grid(igrid)
    use mod_usr_methods, only: usr_var_for_errest, usr_refine_threshold
    use mod_forest, only: coarsen, refine
    use mod_physics, only: phys_energy
    use mod_global_parameters

    integer, intent(in) :: igrid

    integer                            :: iflag, idims, idims2, level
    integer                            :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3, hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3, jxmin1,jxmin2,jxmin3,&
       jxmax1,jxmax2,jxmax3, h2xmin1,h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3,&
        j2xmin1,j2xmin2,j2xmin3,j2xmax1,j2xmax2,j2xmax3, ix1,ix2,ix3
    double precision                   :: epsilon, threshold, wtol(1:nw),&
        xtol(1:ndim)
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3) :: numerator, denominator, error
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: tmp, tmp1, tmp2
    double precision                   :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:nw)
    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)          :: refineflag, coarsenflag

    epsilon = 1.0d-6
    level   = node(plevel_,igrid)
    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
    ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;

    error=zero

    w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)

    if(B0field) then
      if(phys_energy) w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         iw_e)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
         iw_e)+0.5d0*sum(ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,:,0)**2,dim=ndim+1) + sum(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,iw_mag(:))*ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
         ixGlo3:ixGhi3,:,0),dim=ndim+1)
      w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_mag(:))=w(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_mag(:))+ps(igrid)%B0(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,ixGlo3:ixGhi3,:,0)
    end if

    do iflag=1,nw+1

       if(w_refine_weight(iflag)==0.d0) cycle
       numerator=zero

       if (iflag > nw) then
          if (.not. associated(usr_var_for_errest)) then
             call mpistop("usr_var_for_errest not defined")
          else
             call usr_var_for_errest(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
                ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,iflag,ps(igrid)%w,&
                ps(igrid)%x,tmp1)
          end if
       end if

       do idims=1,ndim
          hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
          hxmin3=ixmin3-kr(3,idims);hxmax1=ixmax1-kr(1,idims)
          hxmax2=ixmax2-kr(2,idims);hxmax3=ixmax3-kr(3,idims);
          jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
          jxmin3=ixmin3+kr(3,idims);jxmax1=ixmax1+kr(1,idims)
          jxmax2=ixmax2+kr(2,idims);jxmax3=ixmax3+kr(3,idims);
          if (iflag<=nw) then
            if (logflag(iflag)) then
              tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dlog10(w(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3,iflag))-dlog10(w(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3,iflag))
            else
              tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=w(jxmin1:jxmax1,&
                 jxmin2:jxmax2,jxmin3:jxmax3,iflag)-w(hxmin1:hxmax1,&
                 hxmin2:hxmax2,hxmin3:hxmax3,iflag)
            end if
          else
            if (logflag(iflag)) then
              tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dlog10(tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3))
            else
              tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3)
            end if
          end if
          do idims2=1,ndim
             h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
             h2xmin3=ixMlo3-kr(3,idims2);h2xmax1=ixMhi1-kr(1,idims2)
             h2xmax2=ixMhi2-kr(2,idims2);h2xmax3=ixMhi3-kr(3,idims2);
             j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
             j2xmin3=ixMlo3+kr(3,idims2);j2xmax1=ixMhi1+kr(1,idims2)
             j2xmax2=ixMhi2+kr(2,idims2);j2xmax3=ixMhi3+kr(3,idims2);
             numerator=numerator+(tmp(j2xmin1:j2xmax1,j2xmin2:j2xmax2,&
                j2xmin3:j2xmax3)-tmp(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
                h2xmin3:h2xmax3))**2
          end do
       end do
       denominator=zero
       do idims=1,ndim
          if (iflag<=nw) then
             if (logflag(iflag)) then
              tmp=dabs(dlog10(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
                 iflag)))
             else
              tmp=dabs(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iflag))
             end if
          else
             if (logflag(iflag)) then
              tmp=dabs(dlog10(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
                 ixGlo3:ixGhi3)))
             else
              tmp=dabs(tmp1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3))
             end if
          end if
          hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
          hxmin3=ixmin3-kr(3,idims);hxmax1=ixmax1-kr(1,idims)
          hxmax2=ixmax2-kr(2,idims);hxmax3=ixmax3-kr(3,idims);
          jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
          jxmin3=ixmin3+kr(3,idims);jxmax1=ixmax1+kr(1,idims)
          jxmax2=ixmax2+kr(2,idims);jxmax3=ixmax3+kr(3,idims);
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp(jxmin1:jxmax1,&
             jxmin2:jxmax2,jxmin3:jxmax3)+tmp(hxmin1:hxmax1,hxmin2:hxmax2,&
             hxmin3:hxmax3)
          hxmin1=ixMlo1-2*kr(1,idims);hxmin2=ixMlo2-2*kr(2,idims)
          hxmin3=ixMlo3-2*kr(3,idims);hxmax1=ixMhi1-2*kr(1,idims)
          hxmax2=ixMhi2-2*kr(2,idims);hxmax3=ixMhi3-2*kr(3,idims);
          jxmin1=ixMlo1+2*kr(1,idims);jxmin2=ixMlo2+2*kr(2,idims)
          jxmin3=ixMlo3+2*kr(3,idims);jxmax1=ixMhi1+2*kr(1,idims)
          jxmax2=ixMhi2+2*kr(2,idims);jxmax3=ixMhi3+2*kr(3,idims);
          if (iflag<=nw) then
            if (logflag(iflag)) then
              tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 ixMlo3:ixMhi3)=dabs(dlog10(w(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3,iflag))-dlog10(w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 ixMlo3:ixMhi3,iflag))) +dabs(dlog10(w(ixMlo1:ixMhi1,&
                 ixMlo2:ixMhi2,ixMlo3:ixMhi3,iflag))-dlog10(w(hxmin1:hxmax1,&
                 hxmin2:hxmax2,hxmin3:hxmax3,iflag)))
            else
               tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                  ixMlo3:ixMhi3)=dabs(w(jxmin1:jxmax1,jxmin2:jxmax2,&
                  jxmin3:jxmax3,iflag)-w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                  ixMlo3:ixMhi3,iflag)) +dabs(w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                  ixMlo3:ixMhi3,iflag)-w(hxmin1:hxmax1,hxmin2:hxmax2,&
                  hxmin3:hxmax3,iflag))
            end if
          else
            if (logflag(iflag)) then
              tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 ixMlo3:ixMhi3)=dabs(dlog10(tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3))-dlog10(tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 ixMlo3:ixMhi3))) +dabs(dlog10(tmp1(ixMlo1:ixMhi1,&
                 ixMlo2:ixMhi2,ixMlo3:ixMhi3))-dlog10(tmp1(hxmin1:hxmax1,&
                 hxmin2:hxmax2,hxmin3:hxmax3)))
            else
               tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                  ixMlo3:ixMhi3)=dabs(tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                  jxmin3:jxmax3)-tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                  ixMlo3:ixMhi3)) +dabs(tmp1(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                  ixMlo3:ixMhi3)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
                  hxmin3:hxmax3))
            end if
          end if
          do idims2=1,ndim
             h2xmin1=ixMlo1-kr(1,idims2);h2xmin2=ixMlo2-kr(2,idims2)
             h2xmin3=ixMlo3-kr(3,idims2);h2xmax1=ixMhi1-kr(1,idims2)
             h2xmax2=ixMhi2-kr(2,idims2);h2xmax3=ixMhi3-kr(3,idims2);
             j2xmin1=ixMlo1+kr(1,idims2);j2xmin2=ixMlo2+kr(2,idims2)
             j2xmin3=ixMlo3+kr(3,idims2);j2xmax1=ixMhi1+kr(1,idims2)
             j2xmax2=ixMhi2+kr(2,idims2);j2xmax3=ixMhi3+kr(3,idims2);
             denominator=denominator +(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                ixMlo3:ixMhi3)+amr_wavefilter(level)*(tmp2(j2xmin1:j2xmax1,&
                j2xmin2:j2xmax2,j2xmin3:j2xmax3)+tmp2(h2xmin1:h2xmax1,&
                h2xmin2:h2xmax2,h2xmin3:h2xmax3)))**2
          end do
       end do
       error=error+w_refine_weight(iflag)*dsqrt(numerator/max(denominator,&
          epsilon))
    end do
    
    refineflag=.false.
    coarsenflag=.false.
    threshold=refine_threshold(level)
    do ix3=ixMlo3,ixMhi3
    do ix2=ixMlo2,ixMhi2
    do ix1=ixMlo1,ixMhi1
    
       if (associated(usr_refine_threshold)) then
          wtol(1:nw)   = w(ix1,ix2,ix3,1:nw)
          xtol(1:ndim) = ps(igrid)%x(ix1,ix2,ix3,1:ndim)
          call usr_refine_threshold(wtol, xtol, threshold, global_time, level)
       end if
    
       if (error(ix1,ix2,ix3) >= threshold) then
          refineflag(ix1,ix2,ix3) = .true.
       else if (error(ix1,ix2,ix3) <= derefine_ratio(level)*threshold) then
          coarsenflag(ix1,ix2,ix3) = .true.
       end if
    end do
    end do
    end do

    if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3)).and.level<refine_max_level) refine(igrid,mype)=.true.
    if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3)).and.level>1) coarsen(igrid,mype)=.true.

  end subroutine lohner_grid

  subroutine lohner_orig_grid(igrid)
    use mod_usr_methods, only: usr_var_for_errest, usr_refine_threshold
    use mod_forest, only: coarsen, refine
    use mod_global_parameters

    integer, intent(in) :: igrid

    integer :: iflag, idims, level
    integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxmin1,hxmin2,hxmin3,&
       hxmax1,hxmax2,hxmax3, jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,jxmax3, ix1,&
       ix2,ix3
    double precision :: epsilon, threshold, wtol(1:nw), xtol(1:ndim)
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3) :: numerator, denominator, error
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: dp, dm, dref, tmp1
    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: refineflag, coarsenflag

    epsilon=1.0d-6
    level=node(plevel_,igrid)
    ixmin1=ixMlo1;ixmin2=ixMlo2;ixmin3=ixMlo3;ixmax1=ixMhi1;ixmax2=ixMhi2
    ixmax3=ixMhi3;

    error=zero
    do iflag=1,nw+1
       if(w_refine_weight(iflag)==0.d0) cycle
       numerator=zero
       denominator=zero

       if (iflag > nw) then
          if (.not. associated(usr_var_for_errest)) then
             call mpistop("usr_var_for_errest not defined")
          else
             call usr_var_for_errest(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
                ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,iflag,ps(igrid)%w,&
                ps(igrid)%x,tmp1)
          end if
       end if

       do idims=1,ndim
          hxmin1=ixmin1-kr(1,idims);hxmin2=ixmin2-kr(2,idims)
          hxmin3=ixmin3-kr(3,idims);hxmax1=ixmax1-kr(1,idims)
          hxmax2=ixmax2-kr(2,idims);hxmax3=ixmax3-kr(3,idims);
          jxmin1=ixmin1+kr(1,idims);jxmin2=ixmin2+kr(2,idims)
          jxmin3=ixmin3+kr(3,idims);jxmax1=ixmax1+kr(1,idims)
          jxmax2=ixmax2+kr(2,idims);jxmax3=ixmax3+kr(3,idims);
          if (iflag<=nw) then
            if (logflag(iflag)) then
              dp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dlog10(ps(igrid)%w(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3,iflag))-dlog10(ps(igrid)%w(ixmin1:ixmax1,&
                 ixmin2:ixmax2,ixmin3:ixmax3,iflag))
              dm(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dlog10(ps(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3,iflag))-dlog10(ps(igrid)%w(hxmin1:hxmax1,&
                 hxmin2:hxmax2,hxmin3:hxmax3,iflag))
              dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 ixMlo3:ixMhi3)=dabs(dlog10(ps(igrid)%w(jxmin1:jxmax1,&
                 jxmin2:jxmax2,jxmin3:jxmax3,&
                 iflag)))+ 2.0d0 * dabs(dlog10(ps(igrid)%w(ixMlo1:ixMhi1,&
                 ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                 iflag))) + dabs(dlog10(ps(igrid)%w(hxmin1:hxmax1,&
                 hxmin2:hxmax2,hxmin3:hxmax3,iflag)))
            else
              dp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=ps(igrid)%w(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3,iflag)-ps(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3,iflag)
              dm(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=ps(igrid)%w(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3,iflag)-ps(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3,iflag)
              dref(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
                 ixMlo3:ixMhi3)=dabs(ps(igrid)%w(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3,iflag))+2.0d0*dabs(ps(igrid)%w(ixMlo1:ixMhi1,&
                 ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
                 iflag)) +dabs(ps(igrid)%w(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3,iflag))
            end if
          else
            if (logflag(iflag)) then
              dp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dlog10(tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3))-dlog10(tmp1(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3))
              dm(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dlog10(tmp1(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3))-dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3))
              dref(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dabs(dlog10(tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3)))+ 2.0d0 * dabs(dlog10(tmp1(ixmin1:ixmax1,&
                 ixmin2:ixmax2,ixmin3:ixmax3))) + &
                 dabs(dlog10(tmp1(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3)))
            else
              dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp1(jxmin1:jxmax1,&
                 jxmin2:jxmax2,jxmin3:jxmax3)-tmp1(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)
              dm(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp1(ixmin1:ixmax1,&
                 ixmin2:ixmax2,ixmin3:ixmax3)-tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3)
              dref(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=dabs(tmp1(jxmin1:jxmax1,jxmin2:jxmax2,&
                 jxmin3:jxmax3))+2.0d0*dabs(tmp1(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)) +dabs(tmp1(hxmin1:hxmax1,hxmin2:hxmax2,&
                 hxmin3:hxmax3))
            end if
          end if

          numerator(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)=numerator+(dp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)-dm(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))**2
          denominator(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)=denominator + (dabs(dp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)) + dabs(dm(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)) + amr_wavefilter(level)*dref(ixMlo1:ixMhi1,&
             ixMlo2:ixMhi2,ixMlo3:ixMhi3))**2
  
       end do
       error=error+w_refine_weight(iflag)*dsqrt(numerator/max(denominator,&
          epsilon))
    end do

    refineflag=.false.
    coarsenflag=.false.

    threshold=refine_threshold(level)
    do ix3=ixMlo3,ixMhi3
    do ix2=ixMlo2,ixMhi2
    do ix1=ixMlo1,ixMhi1

       if (associated(usr_refine_threshold)) then
          wtol(1:nw)   = ps(igrid)%w(ix1,ix2,ix3,1:nw)
          xtol(1:ndim) = ps(igrid)%x(ix1,ix2,ix3,1:ndim)
          call usr_refine_threshold(wtol, xtol, threshold, global_time, level)
       end if

       if (error(ix1,ix2,ix3) >= threshold) then
          refineflag(ix1,ix2,ix3) = .true.
       else if (error(ix1,ix2,ix3) <= derefine_ratio(level)*threshold) then
          coarsenflag(ix1,ix2,ix3) = .true.
       end if
    end do
    end do
    end do

    if (any(refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3)).and.level<refine_max_level) refine(igrid,mype)=.true.
    if (all(coarsenflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3)).and.level>1) coarsen(igrid,mype)=.true.

  end subroutine lohner_orig_grid

  subroutine forcedrefine_grid(igrid,w)
    use mod_usr_methods, only: usr_refine_grid
    use mod_forest, only: coarsen, refine, buffer
    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)

    integer :: level
    integer :: my_refine, my_coarsen
    double precision :: qt
    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: refineflag

    level=node(plevel_,igrid)

    ! initialize to 0
    my_refine   = 0
    my_coarsen  = 0

    if (time_advance) then
       qt=global_time+dt
    else
       qt=global_time
    end if

    if (associated(usr_refine_grid)) then
       call usr_refine_grid(igrid,level,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
          ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,qt,w,ps(igrid)%x,&
           my_refine,my_coarsen)
    end if

    if (my_coarsen==1) then
       if (level>1) then
          refine(igrid,mype)=.false.
          coarsen(igrid,mype)=.true.
       else
          refine(igrid,mype)=.false.
          coarsen(igrid,mype)=.false.
       end if
    end if

    if (my_coarsen==-1)then
       coarsen(igrid,mype)=.false.
    end if

    if (my_refine==1) then
       if (level<refine_max_level) then
          refine(igrid,mype)=.true.
          coarsen(igrid,mype)=.false.
       else
          refine(igrid,mype)=.false.
          coarsen(igrid,mype)=.false.
       end if
    end if

    if (my_refine==-1) then
      refine(igrid,mype)=.false.
    end if
  
    if (nbufferx1/=0.or.nbufferx2/=0.or.nbufferx3/=0) then
       if (refine(igrid,mype) .and. .not.buffer(igrid,mype)) then
          refineflag(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)=.true.
          call refinebuffer(igrid,refineflag)
       end if
    end if
  
  end subroutine forcedrefine_grid
  
  subroutine forcedrefine_grid_io(igrid,w)
    use mod_forest, only: coarsen, refine
    use mod_global_parameters

    integer, intent(in)          :: igrid
    double precision, intent(in) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)

    integer                   :: level, my_levmin, my_levmax
    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: refineflag

    level=node(plevel_,igrid)

    if (level_io > 0) then
       my_levmin = level_io
       my_levmax = level_io
    else
       my_levmin = max(1,level_io_min)
       my_levmax = min(refine_max_level,level_io_max)
    end if

    if (level>my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.true.
    elseif (level<my_levmin) then
      refine(igrid,mype)=.true.
      coarsen(igrid,mype)=.false.
    end if

    if (level==my_levmin .or. level==my_levmax) then
      refine(igrid,mype)=.false.
      coarsen(igrid,mype)=.false.
    end if

    if(refine(igrid,mype).and.level>=refine_max_level)refine(igrid,&
       mype)=.false.
    if(coarsen(igrid,mype).and.level<=1)coarsen(igrid,mype)=.false.

  end subroutine forcedrefine_grid_io

  subroutine refinebuffer(igrid,refineflag)
    use mod_forest, only: refine, buffer
    use mod_global_parameters

    integer, intent(in) :: igrid
    logical, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
        intent(in) :: refineflag

    integer :: ishiftbuf1,ishiftbuf2,ishiftbuf3, i1,i2,i3, ixmin1,ixmin2,&
       ixmin3,ixmax1,ixmax2,ixmax3, ineighbor, ipe_neighbor, level

    ishiftbuf1=ixMhi1-ixMlo1-nbufferx1+1;ishiftbuf2=ixMhi2-ixMlo2-nbufferx2+1
    ishiftbuf3=ixMhi3-ixMlo3-nbufferx3+1;
    do i3=-1,1
    do i2=-1,1
    do i1=-1,1
       ixmin1=max(ixMlo1,ixMlo1+i1*ishiftbuf1)
       ixmin2=max(ixMlo2,ixMlo2+i2*ishiftbuf2)
       ixmin3=max(ixMlo3,ixMlo3+i3*ishiftbuf3);
       ixmax1=min(ixMhi1,ixMhi1+i1*ishiftbuf1)
       ixmax2=min(ixMhi2,ixMhi2+i2*ishiftbuf2)
       ixmax3=min(ixMhi3,ixMhi3+i3*ishiftbuf3);
       if (ixmax1<ixmin1.or.ixmax2<ixmin2.or.ixmax3<ixmin3) cycle
       if (any(refineflag(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))) then
          select case (neighbor_type(i1,i2,i3,igrid))
          case (neighbor_coarse)
             ineighbor=neighbor(1,i1,i2,i3,igrid)
             ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
             if (.not.refine(ineighbor,ipe_neighbor)) then
                buffer(ineighbor,ipe_neighbor)=.true.
                refine(ineighbor,ipe_neighbor)=.true.
             end if
          case (neighbor_sibling)
             level=node(plevel_,igrid)
             if (level<refine_max_level) then
                ineighbor=neighbor(1,i1,i2,i3,igrid)
                ipe_neighbor=neighbor(2,i1,i2,i3,igrid)
                if (.not.refine(ineighbor,ipe_neighbor)) then
                   buffer(ineighbor,ipe_neighbor)=.true.
                   refine(ineighbor,ipe_neighbor)=.true.
                end if
             end if
          end select
       end if
    end do
    end do
    end do

  end subroutine refinebuffer

end module mod_errest
