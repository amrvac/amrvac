!> Handles computations for coordinates and variables in output
module mod_calculate_xw
  implicit none

contains

  !> Compute both corner as well as cell-centered values for output
  subroutine calc_grid(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,&
     normconv,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,&
     ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3,first)
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
    use mod_physics, only: phys_energy, physics_type, phys_to_primitive
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: qunit, igrid
    double precision, intent(in), dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,ndim) :: xC
    double precision, intent(in), dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3,ndim)   :: xCC
    integer :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCCmin1,&
       ixCCmin2,ixCCmin3,ixCCmax1,ixCCmax2,ixCCmax3
    logical, intent(in) :: first

    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,ndim) :: xC_TMP
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       ndim)   :: xCC_TMP
    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,nw+nwauxio)   :: wC_TMP
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       nw+nwauxio)     :: wCC_TMP
    double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv

    double precision :: ldw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
        dwC(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision, dimension(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,nw+nwauxio)   :: wC
    double precision, dimension(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
       nw+nwauxio)     :: wCC
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:nw+nwauxio)   :: w
    integer :: nxCC1,nxCC2,nxCC3,idims,jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,&
       jxCmax3,iwe
    integer :: nx1,nx2,nx3, nxC1,nxC2,nxC3, ix1,ix2,ix3, ix, iw, level, idir
    logical, save :: subfirst=.true.

    ! initialize w
    w=0.d0

    ixCmin1=ixMlo1-1;ixCmin2=ixMlo2-1;ixCmin3=ixMlo3-1; ixCmax1=ixMhi1
    ixCmax2=ixMhi2;ixCmax3=ixMhi3; !Corner indices
    ixCCmin1=ixMlo1;ixCCmin2=ixMlo2;ixCCmin3=ixMlo3; ixCCmax1=ixMhi1
    ixCCmax2=ixMhi2;ixCCmax3=ixMhi3; !Center indices

    nx1=ixMhi1-ixMlo1+1;nx2=ixMhi2-ixMlo2+1;nx3=ixMhi3-ixMlo3+1;
    level=node(plevel_,igrid)

    normconv(0) = length_convert_factor
    normconv(1:nw) = w_convert_factor
    block=>ps(igrid)
    w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:nw)=ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw)

    if (nwextra>0) then
     ! here we actually fill the ghost layers for the nwextra variables using
     ! continuous extrapolation (as these values do not exist normally in ghost cells)
     do idims=1,ndim
      select case(idims)
       case(1)
         jxCmin1=ixGhi1+1-nghostcells;jxCmin2=ixGlo2;jxCmin3=ixGlo3;
         jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGhi3;
         do ix1=jxCmin1,jxCmax1
             w(ix1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw) = w(jxCmin1-1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw)
         end do
         jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGlo3;
         jxCmax1=ixGlo1-1+nghostcells;jxCmax2=ixGhi2;jxCmax3=ixGhi3;
         do ix1=jxCmin1,jxCmax1
             w(ix1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw) = w(jxCmax1+1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw)
         end do 
       case(2)
         jxCmin1=ixGlo1;jxCmin2=ixGhi2+1-nghostcells;jxCmin3=ixGlo3;
         jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGhi3;
         do ix2=jxCmin2,jxCmax2
             w(jxCmin1:jxCmax1,ix2,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw) = w(jxCmin1:jxCmax1,jxCmin2-1,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw)
         end do
         jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGlo3;
         jxCmax1=ixGhi1;jxCmax2=ixGlo2-1+nghostcells;jxCmax3=ixGhi3;
         do ix2=jxCmin2,jxCmax2
             w(jxCmin1:jxCmax1,ix2,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw) = w(jxCmin1:jxCmax1,jxCmax2+1,jxCmin3:jxCmax3,&
                nw-nwextra+1:nw)
         end do 
       case(3)
         jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGhi3+1-nghostcells;
         jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGhi3;
         do ix3=jxCmin3,jxCmax3
             w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,ix3,&
                nw-nwextra+1:nw) = w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3-1,&
                nw-nwextra+1:nw)
         end do
         jxCmin1=ixGlo1;jxCmin2=ixGlo2;jxCmin3=ixGlo3;
         jxCmax1=ixGhi1;jxCmax2=ixGhi2;jxCmax3=ixGlo3-1+nghostcells;
         do ix3=jxCmin3,jxCmax3
             w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,ix3,&
                nw-nwextra+1:nw) = w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmax3+1,&
                nw-nwextra+1:nw)
         end do 
      end select
     end do
    end if

    ! next lines needed when usr_aux_output uses gradients
    ! and later on when dwlimiter2 is used
    dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
    dxlevel(3)=rnode(rpdx3_,igrid);
    if(nwauxio>0)then
      ! auxiliary io variables can be computed and added by user
      ! next few lines ensure correct usage of routines like divvector etc
      ! default (no) normalization for auxiliary variables
      normconv(nw+1:nw+nwauxio)=one

      if (.not. associated(usr_aux_output)) then
         call mpistop("usr_aux_output not defined")
      else
         call usr_aux_output(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
            ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,w,&
            ps(igrid)%x,normconv)
      end if
    endif

    ! In case primitives to be saved: use primitive subroutine
    !  extra layer around mesh only needed when storing corner values and averaging
    if(saveprim.and.first) call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
       ixGhi2,ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1,&
       w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nw),ps(igrid)%x)

    if(B0field) then
    ! B0+B1 split handled here
      if(.not.saveprim.and.phys_energy) then
        w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_e)=w(ixGlo1:ixGhi1,&
           ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_e)+&
           0.5d0*sum(ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:,&
           0)**2,dim=ndim+1) + sum(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
           iw_mag(:))*ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:,&
           0),dim=ndim+1)
      end if
      w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_mag(:))=w(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,ixGlo3:ixGhi3,iw_mag(:))+ps(igrid)%B0(ixGlo1:ixGhi1,&
         ixGlo2:ixGhi2,ixGlo3:ixGhi3,:,0)
    end if
    ! compute the cell-center values for w first
    ! cell center values obtained from mere copy
    wCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
       :)=w(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,:)

    ! compute the corner values for w now by averaging

    if(slab_uniform) then
       ! for slab symmetry: no geometrical info required
       do iw=1,nw+nwauxio
         do ix3=ixCmin3,ixCmax3
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
            wC(ix1,ix2,ix3,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,ix3:ix3+1,&
               iw))/dble(2**ndim)
         end do
         end do
         end do
       end do
    else
       do iw=1,nw+nwauxio
         do ix3=ixCmin3,ixCmax3
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
           wC(ix1,ix2,ix3,iw)=sum(w(ix1:ix1+1,ix2:ix2+1,ix3:ix3+1,&
              iw)*ps(igrid)%dvolume(ix1:ix1+1,ix2:ix2+1,&
              ix3:ix3+1)) /sum(ps(igrid)%dvolume(ix1:ix1+1,ix2:ix2+1,&
              ix3:ix3+1))
         end do
         end do
         end do
       end do
    endif

    if(nocartesian) then
      ! keep the coordinate and vector components
      xC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)          = xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)
      wC_TMP(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:nw+nwauxio)    = wC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:nw+nwauxio)
      xCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
         1:ndim)        = xCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
         ixCCmin3:ixCCmax3,1:ndim)
      wCC_TMP(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,ixCCmin3:ixCCmax3,&
         1:nw+nwauxio)  = wCC(ixCCmin1:ixCCmax1,ixCCmin2:ixCCmax2,&
         ixCCmin3:ixCCmax3,1:nw+nwauxio)
    else
      ! do all conversions to cartesian coordinates and vector components
      ! start for the corner values
      call to_cartesian(xC_TMP,wC_TMP,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
         ixCmax3,xC,wC)
      ! then cell center values
      call to_cartesian(xCC_TMP,wCC_TMP,ixCCmin1,ixCCmin2,ixCCmin3,ixCCmax1,&
         ixCCmax2,ixCCmax3,xCC,wCC)
    endif

    ! Warning: differentiate between idl/idlCC/tecplot/tecplotCC/vtu(B)/vtu(B)CC
    if(nwaux>0 .and. mype==0 .and. first.and.subfirst) then
      ! when corner values are computed and auxiliaries present: warn!
      if(convert_type=='idl'.or.convert_type=='tecplot' &
         .or.convert_type=='vtu'.or.convert_type=='vtuB') write(*,&
         *) 'Warning: also averaged auxiliaries within calc_grid'
      subfirst=.false.
    endif

  end subroutine calc_grid

  !> convert to cartesian coordinates and vector components
  subroutine to_cartesian(x_TMP,w_TMP,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
     ixmax3,xC,wC)
    ! conversion of coordinate and vector components from cylindrical/spherical
    ! to cartesian coordinates and components done here
    ! Also: nullifying values lower than smalldouble
    use mod_global_parameters
    use mod_geometry

    integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ix1,ix2,ix3, idim,&
        iw, ivector, iw0
    integer, dimension(nw) :: vectoriw
    double precision :: x_TEC(ndim), w_TEC(nw+nwauxio)
    double precision, dimension(ndim,ndim) :: normal

    double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       ndim) :: xC
    double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       nw+nwauxio)   :: wC

    double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       ndim) :: x_TMP
    double precision, dimension(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       nw+nwauxio)   :: w_TMP

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
    do ix3=ixmin3,ixmax3
    do ix2=ixmin2,ixmax2
    do ix1=ixmin1,ixmax1
       select case (coordinate)
       case (Cartesian,Cartesian_stretched,Cartesian_expansion)
          x_TEC(1:ndim)=xC(ix1,ix2,ix3,1:ndim)
          w_TEC(1:nw+nwauxio)=wC(ix1,ix2,ix3,1:nw+nwauxio)
       case (cylindrical)
          
          
          
          x_TEC(1)=xC(ix1,ix2,ix3,1)*cos(xC(ix1,ix2,ix3,phi_))
          x_TEC(2)=xC(ix1,ix2,ix3,1)*sin(xC(ix1,ix2,ix3,phi_))
          x_TEC(3)=xC(ix1,ix2,ix3,z_)

          if (nvector>0) then
             

             

             
             normal(1,1)=cos(xC(ix1,ix2,ix3,phi_))
             normal(1,phi_)=-sin(xC(ix1,ix2,ix3,phi_))
             normal(1,z_)=zero

             normal(2,1)=sin(xC(ix1,ix2,ix3,phi_))
             normal(2,phi_)=cos(xC(ix1,ix2,ix3,phi_))
             normal(2,z_)=zero

             normal(3,1)=zero
             normal(3,phi_)=zero
             normal(3,z_)=one
          end if
          do iw=1,nw+nwauxio
             if (iw<=nw) iw0=vectoriw(iw)
             if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
                idim=iw-iw0
                w_TEC(iw0+idim)=wC(ix1,ix2,ix3,iw0+1)*normal(idim,1)+wC(ix1,&
                   ix2,ix3,iw0+2)*normal(idim,2)+wC(ix1,ix2,ix3,&
                   iw0+3)*normal(idim,3)
             else
                w_TEC(iw)=wC(ix1,ix2,ix3,iw)
             end if
          end do
       case (spherical)
          x_TEC(1)=xC(ix1,ix2,ix3,1)*sin(xC(ix1,ix2,ix3,2))*cos(xC(ix1,ix2,ix3,&
             3))
          
          
          x_TEC(2)=xC(ix1,ix2,ix3,1)*sin(xC(ix1,ix2,ix3,2))*sin(xC(ix1,ix2,ix3,&
             3))
          x_TEC(3)=xC(ix1,ix2,ix3,1)*cos(xC(ix1,ix2,ix3,2))

          if (nvector>0) then
             
             
             normal(1,1)=sin(xC(ix1,ix2,ix3,2))*cos(xC(ix1,ix2,ix3,3))
             normal(1,2)=cos(xC(ix1,ix2,ix3,2))*cos(xC(ix1,ix2,ix3,3))
             normal(1,3)=-sin(xC(ix1,ix2,ix3,3))

             
             
             normal(2,1)=sin(xC(ix1,ix2,ix3,2))*sin(xC(ix1,ix2,ix3,3))
             normal(2,2)=cos(xC(ix1,ix2,ix3,2))*sin(xC(ix1,ix2,ix3,3))
             normal(2,3)=cos(xC(ix1,ix2,ix3,3))

             normal(3,1)=cos(xC(ix1,ix2,ix3,2))
             normal(3,2)=-sin(xC(ix1,ix2,ix3,2))
             normal(3,3)=zero
          end if
          do iw=1,nw+nwauxio
             if (iw<=nw) iw0=vectoriw(iw)
             if (iw0>=0.and.iw<=iw0+ndim.and.iw<=nw) then
                idim=iw-iw0
                w_TEC(iw0+idim)=wC(ix1,ix2,ix3,iw0+1)*normal(idim,1)+wC(ix1,&
                   ix2,ix3,iw0+2)*normal(idim,2)+wC(ix1,ix2,ix3,&
                   iw0+3)*normal(idim,3)
             else
                w_TEC(iw)=wC(ix1,ix2,ix3,iw)
             end if
          end do
       case default
          write(*,*) "No converter for coordinate=",coordinate
       end select
       x_TMP(ix1,ix2,ix3,1:ndim)=x_TEC(1:ndim)
       w_TMP(ix1,ix2,ix3,1:nw+nwauxio)=w_TEC(1:nw+nwauxio)
       ! Be aware that small values are nullified here!!!
       where(dabs(w_TMP(ix1,ix2,ix3,1:nw+nwauxio))<smalldouble)
             w_TMP(ix1,ix2,ix3,1:nw+nwauxio)=zero
       endwhere
    end do
    end do
    end do

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
    use mod_comm_lib, only: mpistop

    character(len=name_len)   :: wnamei(1:nw+nwauxio),&
       xandwnamei(1:ndim+nw+nwauxio)
    character(len=1024) :: outfilehead

    integer::  space_position,iw,ind
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
          xandwnamei(1)="r"; xandwnamei(2)="Theta"; xandwnamei(3)="Phi"
       case( cylindrical )
          xandwnamei(1)="R";
          
          if( phi_ == 2 )then
             xandwnamei(2)="Phi"
          else
             xandwnamei(2)="Z"
          endif
          
          if( phi_ == 2 )then
             xandwnamei(3)="Z"
          else
             xandwnamei(3)="Phi"
          endif
       case default
          xandwnamei(1)="X"; xandwnamei(2)="Y"; xandwnamei(3)="Z"
    end select

    xandwnamei(ndim+1:ndim+nw+nwauxio)=wnamei(1:nw+nwauxio)

    ! in outfilehead, collect the dimensional names, and all output variable names
    ! first all dimensions
    write(outfilehead,'(a)') TRIM(xandwnamei(1))
    
    do iw=2,ndim
       wname=xandwnamei(iw)
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(&
       wname)
    enddo
   
    ! then all nw variables, with w_write control for inclusion
    do iw=ndim+1,ndim+nw
       wname=xandwnamei(iw)
       if(w_write(iw-ndim)) then
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(&
       wname)
       endif
    enddo
    ! then all nwauxio variables
    if(nwauxio>0) then
      do iw=ndim+nw+1,ndim+nw+nwauxio
         wname=xandwnamei(iw)
    write(outfilehead,'(a)')outfilehead(1:len_trim(outfilehead))//" "//TRIM(&
       wname)
      enddo
    endif

    if(first.and.mype==0)then
      print*,'-------------------------------------------------------------------------------'
      write(unitterm,*&
         )'Saving visual data. Coordinate directions and variable names are:'
      ind=0
      do iw=1,ndim
        ind=ind+1
        print *,ind,xandwnamei(iw)
      enddo
      do iw=1+ndim,nw+ndim
        if(w_write(iw-ndim)) then
          ind=ind+1
          write(*,*) ind,wnamei(iw-ndim)
        end if
      end do
      do iw=ndim+nw+1,ndim+nw+nwauxio
        ind=ind+1
        print *,ind,wnamei(iw-ndim)
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
    double precision, intent(out)     :: xC(ixMlo1-1:ixMhi1,ixMlo2-1:ixMhi2,&
       ixMlo3-1:ixMhi3,ndim)
    double precision, intent(out)     :: xCC(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
       ixMlo3:ixMhi3,ndim)

    integer                           :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
       ixCmax2,ixCmax3, idims, level, ix

    level=node(plevel_,igrid)

    ! coordinates of cell centers
    xCC(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,:)=ps(igrid)%x(ixMlo1:ixMhi1,&
       ixMlo2:ixMhi2,ixMlo3:ixMhi3,:)

    ! coordinates of cell corners
    ixCmin1=ixMlo1-1;ixCmin2=ixMlo2-1;ixCmin3=ixMlo3-1; ixCmax1=ixMhi1
    ixCmax2=ixMhi2;ixCmax3=ixMhi3;
    if(slab_uniform)then
       do idims=1,ndim
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            idims)=ps(igrid)%x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
            idims)+0.5d0*dx(idims,level)
       end do
    else
       ! for any non-cartesian or stretched coordinate (allow multiple stretched directions)
       do ix=ixCmin1,ixCmax1
         xC(ix,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=ps(igrid)%x(ix,&
            ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)+0.5d0*ps(igrid)%dx(ix,&
            ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)
       end do
       do ix=ixCmin2,ixCmax2
         xC(ixCmin1:ixCmax1,ix,ixCmin3:ixCmax3,2)=ps(igrid)%x(ixCmin1:ixCmax1,&
            ix,ixCmin3:ixCmax3,2)+0.5d0*ps(igrid)%dx(ixCmin1:ixCmax1,ix,&
            ixCmin3:ixCmax3,2)
       end do
       do ix=ixCmin3,ixCmax3
         xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ix,3)=ps(igrid)%x(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,ix,3)+0.5d0*ps(igrid)%dx(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,ix,3)
       end do
    endif

  end subroutine calc_x

end module mod_calculate_xw
