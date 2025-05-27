module mod_b0

  implicit none
  private


  public :: set_B0_grid

contains


  subroutine set_B0_grid(igrid)
    use mod_global_parameters
  
    integer, intent(in) :: igrid
  
    integer :: ixCoGmin1,ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3
  
    ixCoGmin1=1;ixCoGmin2=1;ixCoGmin3=1;
    ixCoGmax1=(ixGhi1-2*nghostcells)/2+2*nghostcells
    ixCoGmax2=(ixGhi2-2*nghostcells)/2+2*nghostcells
    ixCoGmax3=(ixGhi3-2*nghostcells)/2+2*nghostcells;
  
    call set_B0_cell(ps(igrid)%B0(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,:,&
       0),ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,ixGlo2,&
       ixGlo3,ixGhi1,ixGhi2,ixGhi3)
    if(B0fieldAllocCoarse) call set_B0_cell(psc(igrid)%B0(ixCoGmin1:ixCoGmax1,&
       ixCoGmin2:ixCoGmax2,ixCoGmin3:ixCoGmax3,:,0),psc(igrid)%x,ixCoGmin1,&
       ixCoGmin2,ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3,ixCoGmin1,ixCoGmin2,&
       ixCoGmin3,ixCoGmax1,ixCoGmax2,ixCoGmax3)
    call set_J0_cell(igrid,ps(igrid)%J0,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
       ixGhi3,ixMlo1-1,ixMlo2-1,ixMlo3-1,ixMhi1+1,ixMhi2+1,ixMhi3+1)
    call set_B0_face(igrid,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
       ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3)
  
  end subroutine set_B0_grid
  
  subroutine set_B0_cell(wB0,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)
    use mod_usr_methods, only: usr_set_B0
    use mod_global_parameters
    use mod_geometry
    
    integer, intent(in):: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
  
    wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1:ndir)=zero
  
    ! approximate cell-averaged B0 as cell-centered B0
    select case (coordinate)
    case (spherical)
       
       if (dabs(Bdip)>smalldouble) then
          wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             1)=2.0d0*Bdip*dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             2))/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**3
          wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             2)=Bdip*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             2))/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**3
       end if
    
       if (abs(Bquad)>smalldouble) then
          wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=wB0(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             1) +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,2)))/x(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,1)**4
          wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)=wB0(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,2)+Bquad*dsin(2.0d0*x(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,1)**4
       end if
       if (abs(Boct)>smalldouble) then
          wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=wB0(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             1) +Boct*(10.0d0*dcos(2.0d0*x(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,2))-2.0d0) *dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             1)**5
          wB0(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)=wB0(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             2) +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,2))) *dsin(x(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,2))/x(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,1)**5
       end if
      
    end select
    
    if (associated(usr_set_B0)) call usr_set_B0(ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,x,&
       wB0)
  
  end subroutine set_B0_cell
  
  subroutine set_J0_cell(igrid,wJ0,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3)
    use mod_usr_methods, only: usr_set_J0
    use mod_global_parameters
    use mod_geometry
  
    integer, intent(in):: igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(inout) :: wJ0(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    integer :: idirmin0, idirmin
  
    if(associated(usr_set_J0)) then
      call usr_set_J0(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixmin1,&
         ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,ps(igrid)%x,wJ0)
    else
      idirmin0 = 7-2*ndir
      call curlvector(ps(igrid)%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,:,0),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,wJ0,idirmin,idirmin0,ndir)
    end if
  
  end subroutine set_J0_cell
  
  subroutine set_B0_face(igrid,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
  
    integer, intent(in) :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
  
    double precision :: delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),xshift1,xshift2,xshift3
    integer :: idims, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ix, idims2
  
    if(slab_uniform)then
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1)=rnode(rpdx1_,&
         igrid)
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2)=rnode(rpdx2_,&
         igrid)
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)=rnode(rpdx3_,&
         igrid)
    else
      ! for all non-cartesian and stretched cartesian coordinates
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=ps(igrid)%dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)
    endif
  
  
    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
      hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
      if(stagger_grid) then
        ! ct needs all transverse cells
        ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
        ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
        ixCmax3=ixOmax3+nghostcells-nghostcells*kr(idims,3)
        ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
        ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2)
        ixCmin3=hxOmin3-nghostcells+nghostcells*kr(idims,3);
      else
        ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
        ixCmin2=hxOmin2;ixCmin3=hxOmin3;
      end if
      ! always xshift=0 or 1/2
      xshift1=half*(one-kr(1,idims));xshift2=half*(one-kr(2,idims))
      xshift3=half*(one-kr(3,idims));
      do idims2=1,ndim
        select case(idims2)
        case(1)
          do ix = ixCmin1,ixCmax1
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ix,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=x(ix,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,1)+(half-xshift1)*delx(ix,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,1)
          end do
        case(2)
          do ix = ixCmin2,ixCmax2
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ixCmin1:ixCmax1,ix,ixCmin3:ixCmax3,2)=x(ixCmin1:ixCmax1,ix,&
               ixCmin3:ixCmax3,2)+(half-xshift2)*delx(ixCmin1:ixCmax1,ix,&
               ixCmin3:ixCmax3,2)
          end do
        case(3)
          do ix = ixCmin3,ixCmax3
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ix,3)=x(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ix,3)+(half-xshift3)*delx(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ix,3)
          end do
        end select
      end do
      call set_B0_cell(ps(igrid)%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,:,idims),xC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3)
    end do
  
  end subroutine set_B0_face


end module mod_b0

