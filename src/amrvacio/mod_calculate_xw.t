!> Handles computations for coordinates and variables in output
module mod_calculate_xw
  implicit none

contains

  !> Compute both corner as well as cell-centered values for output
  subroutine calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                       ixC^L,ixCC^L,first)
    ! this subroutine computes both corner as well as cell-centered values
    ! it handles how we do the center to corner averaging, as well as
    ! whether we switch to cartesian or want primitive or conservative output,
    ! handling the addition of B0 in B0+B1 cases, ...
    !
    ! the normconv is passed on to usr_aux_output for extending with
    ! possible normalization values for the nw+1:nw+nwauxio entries

    use mod_usr_methods, only: usr_aux_output
    use mod_global_parameters
    use mod_limiter
    use mod_physics, only: physics_type, phys_to_primitive

    integer, intent(in) :: qunit, igrid
    double precision, intent(in), dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
    double precision, intent(in), dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC
    integer :: ixC^L,ixCC^L
    logical, intent(in) :: first

    double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
    double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
    double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
    double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
    double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv

    double precision :: ldw(ixG^T), dwC(ixG^T)
    double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC
    double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC
    double precision, dimension(ixG^T,1:nw+nwauxio)   :: w
    integer :: nxCC^D,idims,jxC^L,iwe
    integer :: nx^D, nxC^D, ix^D, ix, iw, level, idir
    logical, save :: subfirst=.true.

    ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D; ! Corner indices
    ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D; ! Center indices

    saveigrid=igrid
    nx^D=ixMhi^D-ixMlo^D+1;
    level=node(plevel_,igrid)

    normconv(0) = length_convert_factor
    normconv(1:nw) = w_convert_factor
    block=>ps(igrid)
    w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)

    if (nwextra>0) then
     ! here we actually fill the ghost layers for the nwextra variables using
     ! continuous extrapolation (as these values do not exist normally in ghost cells)
     do idims=1,ndim
      select case(idims)
       {case(^D)
         jxCmin^DD=ixGhi^D+1-nghostcells^D%jxCmin^DD=ixGlo^DD;
         jxCmax^DD=ixGhi^DD;
         do ix^D=jxCmin^D,jxCmax^D
             w(ix^D^D%jxC^S,nw-nwextra+1:nw) = w(jxCmin^D-1^D%jxC^S,nw-nwextra+1:nw)
         end do
         jxCmin^DD=ixGlo^DD;
         jxCmax^DD=ixGlo^D-1+nghostcells^D%jxCmax^DD=ixGhi^DD;
         do ix^D=jxCmin^D,jxCmax^D
             w(ix^D^D%jxC^S,nw-nwextra+1:nw) = w(jxCmax^D+1^D%jxC^S,nw-nwextra+1:nw)
         end do \}
      end select
     end do
    end if

    ! next lines needed when usr_aux_output uses gradients
    ! and later on when dwlimiter2 is used
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))
    ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
    if(nwauxio>0)then
      ! auxiliary io variables can be computed and added by user
      ! next few lines ensure correct usage of routines like divvector etc
      ! default (no) normalization for auxiliary variables
      normconv(nw+1:nw+nwauxio)=one

      if (.not. associated(usr_aux_output)) then
         call mpistop("usr_aux_output not defined")
      else
         call usr_aux_output(ixG^LL,ixM^LL^LADD1,w,ps(igrid)%x,normconv)
      end if
    endif

    ! In case primitives to be saved: use primitive subroutine
    !  extra layer around mesh only needed when storing corner values and averaging
    if(saveprim.and.first) call phys_to_primitive(ixG^LL,ixM^LL^LADD1,w(ixG^T,1:nw),ps(igrid)%x)

    if(B0field) then
    ! B0+B1 split handled here
      if(.not.saveprim.and.phys_energy) then
        w(ixG^T,iw_e)=w(ixG^T,iw_e)+0.5d0*sum(ps(igrid)%B0(ixG^T,:,0)**2,dim=ndim+1) &
              + sum(w(ixG^T,iw_mag(:))*ps(igrid)%B0(ixG^T,:,0),dim=ndim+1)
      end if
      w(ixG^T,iw_mag(:))=w(ixG^T,iw_mag(:))+ps(igrid)%B0(ixG^T,:,0)
    end if
    ! compute the cell-center values for w first
    ! cell center values obtained from mere copy
    wCC(ixCC^S,:)=w(ixCC^S,:)

    ! compute the corner values for w now by averaging

    if(slab_uniform) then
       ! for slab symmetry: no geometrical info required
       do iw=1,nw+nwauxio
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
            wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw))/dble(2**ndim)
         {end do\}
       end do
    else
       do iw=1,nw+nwauxio
         {do ix^DB=ixCmin^DB,ixCmax^DB\}
           wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw)*ps(igrid)%dvolume(ix^D:ix^D+1)) &
                    /sum(ps(igrid)%dvolume(ix^D:ix^D+1))
         {end do\}
       end do
    endif

    if(nocartesian) then
      ! keep the coordinate and vector components
      xC_TMP(ixC^S,1:ndim)          = xC(ixC^S,1:ndim)
      wC_TMP(ixC^S,1:nw+nwauxio)    = wC(ixC^S,1:nw+nwauxio)
      xCC_TMP(ixCC^S,1:ndim)        = xCC(ixCC^S,1:ndim)
      wCC_TMP(ixCC^S,1:nw+nwauxio)  = wCC(ixCC^S,1:nw+nwauxio)
    else
      ! do all conversions to cartesian coordinates and vector components
      ! start for the corner values
      call to_cartesian(xC_TMP,wC_TMP,ixC^L,xC,wC)
      ! then cell center values
      call to_cartesian(xCC_TMP,wCC_TMP,ixCC^L,xCC,wCC)
    endif

    ! Warning: differentiate between idl/idlCC/tecplot/tecplotCC/vtu(B)/vtu(B)CC
    if(nwaux>0 .and. mype==0 .and. first.and.subfirst) then
      ! when corner values are computed and auxiliaries present: warn!
      if(convert_type=='idl'.or.convert_type=='tecplot' &
       .or.convert_type=='vtu'.or.convert_type=='vtuB') &
          write(*,*) 'Warning: also averaged auxiliaries within calc_grid'
      subfirst=.false.
    endif

  end subroutine calc_grid

  !> convert to cartesian coordinates and vector components
  subroutine to_cartesian(x_TMP,w_TMP,ix^L,xC,wC)
    ! conversion of coordinate and vector components from cylindrical/spherical
    ! to cartesian coordinates and components done here
    ! Also: nullifying values lower than smalldouble
    use mod_global_parameters
    use mod_geometry

    integer :: ix^L, ix^D, idim, iw, ivector, iw0
    integer, dimension(nw) :: vectoriw
    double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)
    double precision, dimension(ndim,ndim) :: normal

    double precision, dimension(ix^S,ndim) :: xC
    double precision, dimension(ix^S,nw+nwauxio)   :: wC

    double precision, dimension(ix^S,ndim) :: x_TMP
    double precision, dimension(ix^S,nw+nwauxio)   :: w_TMP

    iw0=0
    vectoriw=-1
    if(nvector>0) then
      do ivector=1,nvector
         do idim=1,ndim
            vectoriw(iw_vector(ivector)+idim)=iw_vector(ivector)
         end do
      end do
    endif
    x_TEC=0.d0
    {do ix^DB=ixmin^DB,ixmax^DB\}
       select case (coordinate)
       case (Cartesian,Cartesian_stretched,Cartesian_expansion)
          x_TEC(1:ndim)=xC(ix^D,1:ndim)
          w_TEC(1:nw+nwauxio)=wC(ix^D,1:nw+nwauxio)
       case (cylindrical)
          {^IFONED
          x_TEC(1)=xC(ix^D,1)}
          {^IFTWOD
          select case (phi_)
          case (2)
             x_TEC(1)=xC(ix^D,1)*cos(xC(ix^D,2))
             x_TEC(2)=xC(ix^D,1)*sin(xC(ix^D,2))
          case default
             x_TEC(1)=xC(ix^D,1)
             x_TEC(2)=xC(ix^D,2)
          end select}
          {^IFTHREED
          x_TEC(1)=xC(ix^D,1)*cos(xC(ix^D,phi_))
          x_TEC(2)=xC(ix^D,1)*sin(xC(ix^D,phi_))
          x_TEC(3)=xC(ix^D,z_)}

          if (nvector>0) then
             {^IFONED normal(1,1)=one}

             {^IFTWOD
             select case (phi_)
             case (2)
                normal(1,1)=cos(xC(ix^D,2))
                normal(1,2)=-sin(xC(ix^D,2))
                normal(2,1)=sin(xC(ix^D,2))
                normal(2,2)=cos(xC(ix^D,2))
             case default
                normal(1,1)=one
                normal(1,2)=zero
                normal(2,1)=zero
                normal(2,2)=one
             end select}

             {^IFTHREED
             normal(1,1)=cos(xC(ix^D,phi_))
             normal(1,phi_)=-sin(xC(ix^D,phi_))
             normal(1,z_)=zero

             normal(2,1)=sin(xC(ix^D,phi_))
             normal(2,phi_)=cos(xC(ix^D,phi_))
             normal(2,z_)=zero

             normal(3,1)=zero
             normal(3,phi_)=zero
             normal(3,z_)=one}
          end if
          do iw=1,nw+nwauxio
             if (iw<=nw) iw0=vectoriw(iw)
             if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
                idim=iw-iw0
                w_TEC(iw0+idim)={^D&wC(ix^DD,iw0+^D)*normal(idim,^D)+}
             else
                w_TEC(iw)=wC(ix^D,iw)
             end if
          end do
       case (spherical)
          x_TEC(1)=xC(ix^D,1){^NOONED*sin(xC(ix^D,2))}{^IFTHREED*cos(xC(ix^D,3))}
          {^IFTWOD
          x_TEC(2)=xC(ix^D,1)*cos(xC(ix^D,2))}
          {^IFTHREED
          x_TEC(2)=xC(ix^D,1)*sin(xC(ix^D,2))*sin(xC(ix^D,3))
          x_TEC(3)=xC(ix^D,1)*cos(xC(ix^D,2))}

          if (nvector>0) then
             {^IFONED normal(1,1)=one}
             {^NOONED
             normal(1,1)=sin(xC(ix^D,2)){^IFTHREED*cos(xC(ix^D,3))}
             normal(1,2)=cos(xC(ix^D,2)){^IFTHREED*cos(xC(ix^D,3))
             normal(1,3)=-sin(xC(ix^D,3))}}

             {^IFTWOD
             normal(2,1)=cos(xC(ix^D,2))
             normal(2,2)=-sin(xC(ix^D,2))}
             {^IFTHREED
             normal(2,1)=sin(xC(ix^D,2))*sin(xC(ix^D,3))
             normal(2,2)=cos(xC(ix^D,2))*sin(xC(ix^D,3))
             normal(2,3)=cos(xC(ix^D,3))

             normal(3,1)=cos(xC(ix^D,2))
             normal(3,2)=-sin(xC(ix^D,2))
             normal(3,3)=zero}
          end if
          do iw=1,nw+nwauxio
             if (iw<=nw) iw0=vectoriw(iw)
             if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
                idim=iw-iw0
                w_TEC(iw0+idim)={^D&wC(ix^DD,iw0+^D)*normal(idim,^D)+}
             else
                w_TEC(iw)=wC(ix^D,iw)
             end if
          end do
       case default
          write(*,*) "No converter for coordinate=",coordinate
       end select
       x_TMP(ix^D,1:ndim)=x_TEC(1:ndim)
       w_TMP(ix^D,1:nw+nwauxio)=w_TEC(1:nw+nwauxio)
       ! Be aware that small values are nullified here!!!
       where(dabs(w_TMP(ix^D,1:nw+nwauxio))<smalldouble)
             w_TMP(ix^D,1:nw+nwauxio)=zero
       endwhere
    {end do\}

  end subroutine to_cartesian

  !> get all variables names
  subroutine getheadernames(wnamei,xandwnamei,outfilehead)
    ! this collects all variables names in the wnamei character array, getting the info from
    ! the prim_wnames/cons_wnames strings (depending on saveprim). It combines this info with names
    ! for the dimensional directions in the xandwnamei array. In the outfilehead, it collects
    ! the dimensional names, and only those names from the nw variables for output (through w_write)
    ! together with the added names for nwauxio variables
    use mod_usr_methods, only: usr_add_aux_names
    use mod_global_parameters
    use mod_geometry

    character(len=name_len)   :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
    character(len=1024) :: outfilehead

    integer::  space_position,iw
    character(len=name_len)::  wname
    character(len=std_len):: aux_variable_names
    character(len=std_len)::  scanstring

    logical, save:: first=.true.

    ! in case additional variables are computed and stored for output
    if (nwauxio>0) then
       if (.not. associated(usr_add_aux_names)) then
          call mpistop("usr_add_aux_names not defined")
       else
          call usr_add_aux_names(aux_variable_names)
       end if
    end if

    ! --- part to provide variable names from prim_wnames/cons_wnames strings
    if(saveprim) then
       scanstring=TRIM(aux_variable_names)
       wnamei(1:nw)=prim_wnames(1:nw)
    else
       scanstring=TRIM(aux_variable_names)
       wnamei(1:nw)=cons_wnames(1:nw)
    endif

    space_position=index(scanstring,' ')
    do iw=nw+1,nw+nwauxio
       do while (space_position==1)
         scanstring=scanstring(2:)
         space_position=index(scanstring,' ')
       enddo
       wname=scanstring(:space_position-1)
       scanstring=scanstring(space_position+1:)
       space_position=index(scanstring,' ')

       ! fill all names, even those that we will not write away from the first nw
       wnamei(iw)=TRIM(wname)
    enddo
    ! --- end of part to provide variable names

    select case (coordinate)
       case( spherical )
          xandwnamei(1)="r";{^NOONED xandwnamei(2)="Theta"};{^IFTHREED xandwnamei(3)="Phi"}
       case( cylindrical )
          xandwnamei(1)="R";
          {^NOONED
          if( phi_ == 2 )then
             xandwnamei(2)="Phi"
          else
             xandwnamei(2)="Z"
          endif}
          {^IFTHREED
          if( phi_ == 2 )then
             xandwnamei(3)="Z"
          else
             xandwnamei(3)="Phi"
          endif}
       case default
          xandwnamei(1)="X";{^NOONED xandwnamei(2)="Y"};{^IFTHREED xandwnamei(3)="Z"}
    end select

    xandwnamei(ndim+1:ndim+nw+nwauxio)=wnamei(1:nw+nwauxio)

    ! in outfilehead, collect the dimensional names, and all output variable names
    ! first all dimensions
    write(outfilehead,'(a)') TRIM(xandwnamei(1))
    {^NOONED
    do iw=2,ndim
       wname=xandwnamei(iw)
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
    enddo
    }
    ! then all nw variables, with w_write control for inclusion
    do iw=ndim+1,ndim+nw
       wname=xandwnamei(iw)
       if(w_write(iw-ndim)) then
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
       endif
    enddo
    ! then all nwauxio variables
    if(nwauxio>0) then
      do iw=ndim+nw+1,ndim+nw+nwauxio
         wname=xandwnamei(iw)
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(wname)
      enddo
    endif

    if(first.and.mype==0)then
      print*,'-------------------------------------------------------------------------------'
      write(unitterm,*)'Saving visual data. Coordinate directions and variable names are:'
      do iw=1,ndim
        print *,iw,xandwnamei(iw)
      enddo
      do iw=ndim+1,ndim+nw+nwauxio
        print *,iw,wnamei(iw-ndim),xandwnamei(iw)
      enddo
      write(unitterm,*)'time =', global_time
      print*,'-------------------------------------------------------------------------------'
      first=.false.
    endif

  end subroutine getheadernames

  !> computes cell corner (xC) and cell center (xCC) coordinates
  subroutine calc_x(igrid,xC,xCC)
    use mod_global_parameters

    integer, intent(in)               :: igrid
    double precision, intent(out)     :: xC(ixMlo^D-1:ixMhi^D,ndim)
    double precision, intent(out)     :: xCC(ixMlo^D:ixMhi^D,ndim)

    integer                           :: ixC^L, idims, level, ix

    level=node(plevel_,igrid)

    ! coordinates of cell centers
    xCC(ixM^T,:)=ps(igrid)%x(ixM^T,:)

    ! coordinates of cell corners
    ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D;
    if(slab_uniform)then
       do idims=1,ndim
         xC(ixC^S,idims)=ps(igrid)%x(ixC^S,idims)+0.5d0*dx(idims,level)
       end do
    else
       ! for any non-cartesian or stretched coordinate (allow multiple stretched directions)
       {do ix=ixCmin^D,ixCmax^D
         xC(ix^D%ixC^S,^D)=ps(igrid)%x(ix^D%ixC^S,^D)+0.5d0*ps(igrid)%dx(ix^D%ixC^S,^D)
       end do\}
    endif

  end subroutine calc_x

end module mod_calculate_xw
