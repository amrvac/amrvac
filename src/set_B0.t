subroutine set_B0_grid(igrid)
  use mod_global_parameters
  use mod_mhd_phys, only: mhd_semirelativistic

  integer, intent(in) :: igrid

  integer :: ixCoG^L

  ixCoGmin^D=1;
  ixCoGmax^D=(ixGhi^D-2*nghostcells)/2+2*nghostcells;

  call set_B0_cell(ps(igrid)%B0(ixG^T,:,0),ps(igrid)%x,ixG^LL,ixG^LL)
  if(mhd_semirelativistic) call set_B0_cell(psc(igrid)%B0(ixCoG^S,:,0),psc(igrid)%x,ixCoG^L,ixCoG^L)
  call set_J0_cell(igrid,ps(igrid)%J0,ixG^LL,ixM^LL^LADD1)
  call set_B0_face(igrid,ps(igrid)%x,ixG^LL,ixM^LL)

end subroutine set_B0_grid

subroutine set_B0_cell(wB0,x,ixI^L,ix^L)
  use mod_usr_methods, only: usr_set_B0
  use mod_global_parameters
  use mod_geometry
  
  integer, intent(in):: ixI^L,ix^L
  double precision, intent(inout) :: wB0(ixI^S,1:ndir)
  double precision, intent(in) :: x(ixI^S,1:ndim)

  wB0(ix^S,1:ndir)=zero

  ! approximate cell-averaged B0 as cell-centered B0
  select case (coordinate)
  case (spherical)
     {^NOONED
     if (dabs(Bdip)>smalldouble) then
        wB0(ix^S,1)=2.0d0*Bdip*dcos(x(ix^S,2))/x(ix^S,1)**3
        wB0(ix^S,2)=Bdip*dsin(x(ix^S,2))/x(ix^S,1)**3
     end if
  
     if (abs(Bquad)>smalldouble) then
        wB0(ix^S,1)=wB0(ix^S,1) &
             +Bquad*0.5d0*(1.0d0+3.0d0*dcos(2.0d0*x(ix^S,2)))/x(ix^S,1)**4
        wB0(ix^S,2)=wB0(ix^S,2)+Bquad*dsin(2.0d0*x(ix^S,2))/x(ix^S,1)**4
     end if
     if (abs(Boct)>smalldouble) then
        wB0(ix^S,1)=wB0(ix^S,1) &
                     +Boct*(10.0d0*dcos(2.0d0*x(ix^S,2))-2.0d0) &
                          *dcos(x(ix^S,2))/x(ix^S,1)**5
        wB0(ix^S,2)=wB0(ix^S,2) &
                     +Boct*1.5d0*(3.0d0+5.0d0*dcos(2.0d0*x(ix^S,2))) &
                          *dsin(x(ix^S,2))/x(ix^S,1)**5
     end if
     }
  end select
  
  if (associated(usr_set_B0)) call usr_set_B0(ixI^L,ix^L,x,wB0)

end subroutine set_B0_cell

subroutine set_J0_cell(igrid,wJ0,ixI^L,ix^L)
  use mod_usr_methods, only: usr_set_J0
  use mod_global_parameters
  use mod_geometry

  integer, intent(in):: igrid,ixI^L,ix^L
  double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:3)
  integer :: idirmin0, idirmin

  if(associated(usr_set_J0)) then
    call usr_set_J0(ixI^L,ix^L,ps(igrid)%x,wJ0)
  else
    idirmin0 = 7-2*ndir
    call curlvector(ps(igrid)%B0(ixI^S,:,0),ixI^L,ix^L,wJ0,idirmin,idirmin0,ndir)
  end if

end subroutine set_J0_cell

subroutine set_B0_face(igrid,x,ixI^L,ixO^L)
  use mod_global_parameters

  integer, intent(in) :: igrid, ixI^L, ixO^L
  double precision, intent(in) :: x(ixI^S,1:ndim)

  double precision :: delx(ixI^S,1:ndim)
  double precision :: xC(ixI^S,1:ndim),xshift^D
  integer :: idims, ixC^L, hxO^L, ix, idims2

  if(slab_uniform)then
    ^D&delx(ixI^S,^D)=rnode(rpdx^D_,igrid)\
  else
    ! for all non-cartesian and stretched cartesian coordinates
    delx(ixI^S,1:ndim)=ps(igrid)%dx(ixI^S,1:ndim)
  endif


  do idims=1,ndim
    hxO^L=ixO^L-kr(idims,^D);
    if(stagger_grid) then
      ! ct needs all transverse cells
      ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
    else
      ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
      ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
    end if
    ! always xshift=0 or 1/2
    xshift^D=half*(one-kr(^D,idims));
    do idims2=1,ndim
      select case(idims2)
      {case(^D)
        do ix = ixC^LIM^D
          ! xshift=half: this is the cell center coordinate
          ! xshift=0: this is the cell edge i+1/2 coordinate
          xC(ix^D%ixC^S,^D)=x(ix^D%ixC^S,^D)+(half-xshift^D)*delx(ix^D%ixC^S,^D)
        end do\}
      end select
    end do
    call set_B0_cell(ps(igrid)%B0(ixI^S,:,idims),xC,ixI^L,ixC^L)
  end do

end subroutine set_B0_face
